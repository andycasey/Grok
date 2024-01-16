
import numpy as np


import numpy as np
import gzip
import warnings
import pickle
import os
from itertools import cycle
from functools import cache
from typing import Optional, Tuple, Sequence, Union
from scipy import optimize as op
from sklearn.decomposition._nmf import _fit_coordinate_descent
from sklearn.exceptions import ConvergenceWarning
from scipy.signal.windows import tukey
from scipy.optimize import leastsq
from scipy import interpolate


from specutils.lsf import instrument_lsf_kernel

def expand_path(path):
    return os.path.abspath(os.path.expanduser(path))

@cache
def load_basis_vectors(path, P, pad=0):    
    full_path = expand_path(path)
    if full_path.lower().endswith(".gz"):
        with gzip.open(full_path, "rb") as fp:
            masked_basis_vectors = pickle.load(fp)
    else:        
        with open(full_path, "rb") as fp:
            masked_basis_vectors = pickle.load(fp)
    if pad > 0:
        basis_vectors = np.zeros((masked_basis_vectors.shape[0], P + 2 * pad))
        basis_vectors[:, pad:-pad] = masked_basis_vectors
    else:
        basis_vectors = masked_basis_vectors
    return basis_vectors



class ContinuumModel:
    
    """
    Simultaneously model the continuum in multiple spectra of the same star,
    using a pre-trained basis set of spectra, or data-driven.
    """
    
    def __init__(
        self,
        basis_vacuum_wavelength=None,
        basis_vectors=None,
    ):
        self.basis_vacuum_wavelength = basis_vacuum_wavelength
        self.basis_vectors = basis_vectors
        
        
    def fit(
        self,
        wavelength: Sequence[Sequence[float]], 
        flux: Sequence[Sequence[float]], 
        ivar: Sequence[Sequence[float]], 
        deg: int, #Union[int, Sequence[int]],
        L: float, #Union[float, Sequence[float]],
        vacuum: Optional[bool] = True,
        alpha: Optional[Union[float, int]] = 0,
        ax=None
    ):
        C, P = self.basis_vectors.shape
        
        O = len(wavelength)        
        #deg, L = (_cycle(deg), _cycle(L))
        
        F = 2 * deg + 1
        #F = deg + 1
        
        # Interpolate the basis vectors to the observed wavelengths.
        if vacuum:
            basis_wavelength = self.basis_vacuum_wavelength
        else:
            basis_wavelength = vacuum_to_air(self.basis_vacuum_wavelength)
        
        wavelength, flux, ivar = _ensure_ragged(wavelength, flux, ivar)
        
        # Make ragged
        order_indices, wls, fluxs, ivars = ([], [], [], [])
        for o, (wl_, flux_, ivar_) in enumerate(zip(*(wavelength, flux, ivar))):
            order_indices.extend([o] * len(wl_))
            fluxs.extend(flux_)
            ivars.extend(ivar_)
            wls.extend(wl_)
        
        order_indices, fluxs, ivars, wls = tuple(map(np.array, (order_indices, fluxs, ivars, wls)))
        indices = np.argsort(wls)
        wls, fluxs, ivars, order_indices = tuple(map(lambda x: x[indices], (wls, fluxs, ivars, order_indices)))
        
        ivars[ivars < 1e-10] = 0
        
        print("ASSUMING R_init and R_new")
        from time import time
        t_init = time()
        K = instrument_lsf_kernel(basis_wavelength, wls, 3000, 85000)
        t_before = time() - t_init
        print(t_before)
        basis_vectors = self.basis_vectors @ K
        
        A = np.zeros((wls.size, F * O))
        p0 = np.hstack([10**(-2) * np.ones(C), np.zeros(O * F)])
        
        for o in range(O):
            mask = (order_indices == o)
            if L is None:
                L_ = 2 * np.ptp(wavelength[o])
            else:
                L_ = L
            As = design_matrix(wls[mask], deg, L_).T
            #As = np.vander(wls[mask], deg + 1)
            A[mask, o*F:(o+1)*F] = As
            (continuum, theta, _) = fit_sines_and_cosines(wls[mask], fluxs[mask], ivars[mask], As)
            #theta = _fit_polynomial(wls[mask], fluxs[mask], ivars[mask], deg)
            p0[C + o*F:C + (o+1)*F] = theta
        
        def f(_, *theta):
            print(theta)
            W = theta[:C]
            return (1 - W @ basis_vectors) * (A @ theta[C:])
            
        p_opt, cov = op.curve_fit(
            f,
            None,
            fluxs,
            p0=p0,
            sigma=ivars**-0.5,
            bounds=self.get_bounds(p0.size)
        )
        
        y_pred = f(None, *p_opt)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
            
        rectified_model_flux = 1 - p_opt[:C] @ basis_vectors
        continuum = A @ p_opt[C:]
        for i in range(O):
            mask = (order_indices == i)
            #ls = ":" if fig is None else "-"
            ls = "-"
            mean = np.mean(wavelength[i])
            ax.plot(wls - mean, i + rectified_model_flux, c="tab:red", ls=ls)
            ax.plot(wavelength[i] - mean, i + flux[i] / continuum[mask], c='k', ls=ls)
            #sw, ew = wls[mask][[0, -1]]
            #ax.set_xlim(sw - 10, ew + 10)
        return fig
    
    def get_bounds(self, N, component_bounds=(0, +np.inf)):
        C, P = self.basis_vectors.shape          

        return np.vstack([
            np.tile(component_bounds, C).reshape((C, 2)),
            np.tile([-np.inf, +np.inf], N - C).reshape((-1, 2))
        ]).T            

def _fit_polynomial(wavelength, flux, ivar, deg):
    
    A = np.vander(wavelength, deg + 1)
    C_inv = np.diag(ivar)
    MTM = A.T @ C_inv @ A
    MTy = A.T @ C_inv @ flux[:, None]
    return np.linalg.solve(MTM, MTy).flatten()
    
def fit_sines_and_cosines(wavelength, flux, ivar, A = None, **kwargs):
    if A is None:
        A = design_matrix(wavelength, **kwargs)
    MTM = A.T @ (ivar[:, None] * A)
    MTy = A.T @ (ivar * flux)
    try:
        theta = np.linalg.solve(MTM, MTy)
    except np.linalg.LinAlgError:
        if np.any(ivar > 0):
            raise
    continuum = A @ theta
    return (continuum, theta, A)    
        
def design_matrix(wavelength: np.array, deg: int, L) -> np.array:
    #L = 2 * np.ptp(wavelength)
    scale = 2 * (np.pi / L)
    return np.vstack(
        [
            np.ones_like(wavelength).reshape((1, -1)),
            np.array(
                [
                    [np.cos(o * scale * wavelength), np.sin(o * scale * wavelength)]
                    for o in range(1, deg + 1)
                ]
            ).reshape((2 * deg, wavelength.size)),
        ]
    )    
    
        

def _ensure_ragged(wavelength, flux, ivar):
    if isinstance(wavelength[0], (float, int)):
        return ([wavelength], [flux], [ivar])
    return (wavelength, flux, ivar)


def air_to_vacuum(lambdas):
    """
    Convert air wavelengths to vacuum.

    As per: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    
    :param lambdas:
        The air wavelengths.
    """    
    s = 10**4 / lambdas
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    return lambdas * n 

def vacuum_to_air(lambdas):
    """
    Convert vacuum wavelengths to air.

    The formula from Donald Morton (2000, ApJ. Suppl., 130, 403) is used for the 
    refraction index, which is also the IAU standard.

    As per https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

    :param lambdas:
        The vacuum wavelengths.
    """
    s = 10**4 / lambdas
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
    return (lambdas / n) 



if __name__ == "__main__":
    
    
    import gzip
    import pickle
    
    basis_vacuum_wavelength = 10 * (10**(2.57671464 + np.arange(167283) * 2.25855074e-06))

    with gzip.open("../../H.pkl.gz", "rb") as fp:
        basis_vectors = pickle.load(fp)
        
    from specutils import Spectrum1D
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563blue_multi.fits", flux_ext=6)
    import matplotlib.pyplot as plt
    trim_blue = 0
    v_rel = -42.3
    wl, flux, ivar = ([], [], [])
    for index in range(0, len(spectra)):
        wl.append(spectra[index].wavelength[trim_blue:] * (1 - v_rel/299792.458))
        flux.append(spectra[index].flux[trim_blue:])
        ivar.append(spectra[index].ivar[trim_blue:])    
    '''
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563red_multi.fits", flux_ext=6)
    for index in range(0, len(spectra)):
        wl.append(spectra[index].wavelength[trim_blue:] * (1 - v_rel/299792.458))
        flux.append(spectra[index].flux[trim_blue:])
        ivar.append(spectra[index].ivar[trim_blue:])        
    '''
    raise a
    model = ContinuumModel(basis_vacuum_wavelength, basis_vectors)
    
    from astropy.io import fits
    image = fits.open("/Users/andycasey/research/metal-free-stars/spectra/lamost/331645249907140480.fits")
    
    wl = [image[1].data["WAVELENGTH"][0, 1:]]
    flux = [image[1].data["FLUX"][0, 1:]]
    ivar = [image[1].data["IVAR"][0, 1:]]
    
    fig = model.fit(wl, flux, ivar, deg=3, L=10000, vacuum=False)
    #fig.axes[0].set_title(f"index = {index}")
    
    
    raise a
    for index in range(21, 35):
        
        wl = [spectra[index].wavelength[trim_blue:] * (1 - v_rel/299792.458)]
        flux = [spectra[index].flux[trim_blue:]]
        ivar = [spectra[index].ivar[trim_blue:]       ]
        
        model = ContinuumModel(basis_vacuum_wavelength, basis_vectors)
        fig = model.fit(wl, flux, ivar, deg=3, L=None, vacuum=False)
        fig.axes[0].set_title(f"index = {index}")