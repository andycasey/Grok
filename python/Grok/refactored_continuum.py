
import numpy as np

import h5py as h5

import numpy as np
#import autograd.numpy as np
#from autograd import jacobian as auto_jacobian # consider make_jvp
import gzip
import warnings
import pickle
import os
from time import time

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

from threading import Thread

class CallbackThread(Thread):
    def __init__(self, args=(), callback=None, **kwargs):
        target = kwargs.pop("target")
        super(CallbackThread, self).__init__(target=self.target_with_callback, **kwargs)
        self.callback = callback
        self.method = target
        self.args = args
        return None
    
    def target_with_callback(self):
        result = self.method(*self.args)
        if self.callback is not None:
            self.callback(result)

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


def _convolve_basis_vectors(basis_vectors, λi, λo, Ro, Ri):
    return basis_vectors @ instrument_lsf_kernel(λi, λo, Ro, Ri)


def _interpolate_basis_vectors(basis_vectors, λi, λo):        
    bv = np.zeros((basis_vectors.shape[0], len(λo)))
    for c, basis_vector in enumerate(basis_vectors):
        bv[c] = np.interp(λo, λi, basis_vector, left=0, right=0)
    return bv

def get_L_default(wavelength):
    return 2 * np.ptp(wavelength)


class ContinuumBasis:
    pass

class Sinusoids(ContinuumBasis):    
        
    def __init__(self, deg=3, L=None):
        self.deg = deg
        self.L = L
        return None
    
    @property
    def num_parameters(self):
        return 2 * self.deg + 1
        
    def default_L(self, wavelength):
        return 2 * np.ptp(wavelength)
    
    def design_matrix(self, wavelength):
        L = self.L or self.default_L(wavelength)
        scale = 2 * (np.pi / L)
        return np.vstack(
            [
                np.ones_like(wavelength).reshape((1, -1)),
                np.array(
                    [
                        [np.cos(o * scale * wavelength), np.sin(o * scale * wavelength)]
                        for o in range(1, self.deg + 1)
                    ]
                ).reshape((2 * self.deg, wavelength.size)),
            ]
        ).T

class Polynomial(ContinuumBasis):
    
    def __init__(self, deg=2):
        self.deg = deg
        return None
    
    @property
    def num_parameters(self):
        return self.deg + 1
    
    def design_matrix(self, wavelength):
        return np.vander(wavelength, self.deg + 1)


class ContinuumModel:
    
    """Model the continuum in many spectra of the a star using pre-trained basis vectors."""
    
    def __init__(self, basis_vacuum_wavelength, basis_vectors, basis_meta=None):
        if basis_meta is None:
            basis_meta = dict()
        self.basis_vacuum_wavelength = basis_vacuum_wavelength
        self.basis_vectors = basis_vectors            
        self.basis_meta = basis_meta    
        return None        
    
    
    @classmethod
    def from_path(cls, path):
        """
        Load a continuum model from a file.
        
        :param path:
            The path to the continuum model file.
        """
        with h5.File(expand_path(path), "r") as fp:
            basis_vacuum_wavelength = fp["vacuum_wavelength"][:]
            basis_vectors = fp["basis_vectors"][:]
            basis_meta = dict(fp["basis_vectors"].attrs)            
        return cls(
            basis_vacuum_wavelength,
            basis_vectors,
            basis_meta
        )
    
    
    def resample_basis_vectors(self, λ, Ro=None, Ri=None, vacuum=True, asynchronous=False):
        """
        Resample (and optionally convolve) the basis vectors to the observed wavelengths.
        
        :param λ:
            The observed rest-frame wavelengths.
        
        :param Ro: [optional]
            The spectral resolution of the observations. If `None` is given, then no
            convolution will take place; the basis vectors will be interpolated.
        
        :param Ri: [optional]
            The spectral resolution of the basis vectors. If `None` is given then
            it defaults to the input spectral resolution stored in the metadata of
            the basis vector file, or infinity if the input spectral resolution is
            not stored in that file.
        
        :param vacuum: [optional]
            Specify whether the observed rest-frame wavelengths are in vacuum (True) or 
            air (False).
        
        :param asynchronous: [optional]
            If true, the convolution will be performed in a separate thread.
        """
        callback = self._resample_basis_vectors_complete
        f = (lambda x: x) if vacuum else vacuum_to_air        
        if Ro is None:
            target = _interpolate_basis_vectors
            args = (self.basis_vectors, f(self.basis_vacuum_wavelength), λ)
        else:
            target = _convolve_basis_vectors
            Ri = Ri or self.basis_meta.get("Ri", np.inf)
            args = (self.basis_vectors, f(self.basis_vacuum_wavelength), λ, Ro, Ri)

        if asynchronous:
            return CallbackThread(callback=callback, target=target, args=args).start()
        else:
            return callback(target(*args))
        
    def _resample_basis_vectors_complete(self, basis_vectors):
        self.resampled_basis_vectors = basis_vectors
        return basis_vectors


    def get_initial_guess(self, A, flux, ivar, order_masks, N, p0=None, small=1e-10):
        C, _ = self.basis_vectors.shape
        x0 = small * np.ones(C + A.shape[1])
        if p0 is None:
            scalar = np.ones_like(flux)
        else:
            scalar = 1 - p0 @ self.resampled_basis_vectors
            x0[:C] = p0
        si = 0
        for n, mask in zip(N, order_masks):      
            scale = scalar[mask]      
            x0[C + si:C + si + n] = _solve_X(flux[mask] / scale, ivar[mask] * scale**2, A[mask, si:si+n])
            si += n
        return x0
    
    
    def fit(
        self,
        λ: Sequence[Sequence[float]], 
        flux: Sequence[Sequence[float]], 
        ivar: Sequence[Sequence[float]], 
        continuum_basis: Union[ContinuumBasis, Sequence[ContinuumBasis]] = Sinusoids,
        R: Optional[float] = None,
        vacuum: Optional[bool] = True,        
        p0: Optional[Sequence[float]] = None,
        alpha: Optional[Union[float, int]] = 0,
        ivar_min: Optional[float] = 1e-8,
        callback_iteration: Optional[callable] = None,
        callback_complete: Optional[callable] = None
    ):
        """
        Fit the model to some observed spectra.
        
        :param λ:

        """
        C, _ = self.basis_vectors.shape
        
        # Prepare data arrays.
        O = len(λ)
        λ, flux, ivar = _ensure_ragged(O, λ, flux, ivar)
        
        oi = np.hstack([np.ones_like(wl) * o for o, wl in enumerate(λ)]) # order indices
        λ = np.hstack(λ)
        pi = np.argsort(λ) # pixel indices
        
        λ, flux, ivar, oi = (λ[pi], np.hstack(flux)[pi], np.hstack(ivar)[pi], oi[pi])
        
        bad_pixel = (~np.isfinite(flux)) | (~np.isfinite(ivar)) | (ivar <= ivar_min)
        ivar[bad_pixel] = 0
        
        # Prepare basis vectors.
        try:
            # Use pre-convolution if we have it.
            basis_vectors = self.resampled_basis_vectors 
        except AttributeError:
            basis_vectors = self.resample_basis_vectors(λ, Ro=R, vacuum=vacuum)
        else:
            if R is not None:
                warnings.warn("R was given but the basis vectors have already been convolved; ignoring R")
        
        # Compute number of continuum parameters per order.
        continuum_basis = _expand_continuum_basis(continuum_basis, O)
        
        N = [b.num_parameters for b in continuum_basis]
                
        A = np.zeros((λ.size, sum(N)))
        
        # Construct the full design matrix.
        si, order_masks = (0, [])
        for o, (cb, n) in enumerate(zip(continuum_basis, N)):
            mask = (oi == o)
            order_masks.append(mask)            
            A[mask, si:si+n] = cb.design_matrix(λ[mask])
            si += n

        if callback_iteration is None:
            def f(theta):
                return ((1 - theta[:C] @ basis_vectors) * (A @ theta[C:]) - flux) * inv_sigma
        else:
            def f(theta):
                y = (1 - theta[:C] @ basis_vectors) * (A @ theta[C:])
                callback_iteration(theta, y)
                return (y - flux) * inv_sigma
            
        # Get an initial guess of continuum with no absorption.
        x0 = self.get_initial_guess(A, flux, ivar, order_masks, N, p0)
        
        inv_sigma = np.sqrt(ivar)
        
        # Let's pre-compute some things for faster Jacobian evaluations.
        mBVT = -basis_vectors.T
        A_inv_sigma = A * inv_sigma[:, None]
        mBVT_inv_sigma = mBVT * inv_sigma[:, None]
        
        def jacobian(theta):
            return np.hstack([
                (A @ theta[C:, None]) * mBVT_inv_sigma,
                (1 + mBVT @ theta[:C, None]) * A_inv_sigma
            ])
        
        print("Starting")
        t_init = time()
        result = op.least_squares(
            f, x0, jac=jacobian, bounds=self.get_bounds(x0.size),
            verbose=2
        )
        t_op = time() - t_init
        print(f"took {t_op:.1f} s to optimize {len(x0)} parameters with {flux.size} data points")
        
        y_pred = np.array([A @ result.x[C:], 1 - result.x[:C] @ basis_vectors])                
        continuum, rectified_model_flux = zip(*[y_pred[:, m] for m in order_masks])

        return (result.x, None, continuum, rectified_model_flux)
    
    
    def get_bounds(self, N, component_bounds=(1e-10, +np.inf)):
        C, _ = self.basis_vectors.shape          
        return np.vstack([
            np.tile(component_bounds, C).reshape((C, 2)),
            np.tile([-np.inf, +np.inf], N - C).reshape((-1, 2))
        ]).T            


def _solve_X(flux: Sequence[float], ivar: Sequence[float], A: np.array):
    MTM = A.T @ (ivar[:, None] * A)
    MTy = A.T @ (ivar * flux)
    theta = np.linalg.solve(MTM, MTy)
    return theta
        
        
def instantiate(item, **kwargs):
    if isinstance(item, type):
        return item(**kwargs)
    else:
        return item

def _expand_continuum_basis(continuum_basis, N):
    if isinstance(continuum_basis, (list, tuple)):
        if len(continuum_basis) != N:
            raise ValueError(f"a sequence of continuum_basis was given, but the length does not match the number of orders ({len(continuum_basis)} != {N})")
        return tuple(map(instantiate, continuum_basis))
    else:
        return tuple([instantiate(continuum_basis)] * N)

def _ensure_ragged(N, *args):
    if isinstance(args[0], (list, tuple, np.array)):
        if len(args[0]) != N:
            raise ValueError(f"the length of {args[0]} must match {N}")
        return args
    else:
        return tuple(map(list, args))


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
    
    
    model = ContinuumModel.from_path("bosz-highres-optical-basis-vectors.h5")
    
    '''
    from astropy.io import fits
    image = fits.open("/Users/andycasey/research/metal-free-stars/spectra/lamost/331645249907140480.fits")
    
    wl = [image[1].data["WAVELENGTH"][0, 1:] * (1 - 35/299792.458)]
    flux = [image[1].data["FLUX"][0, 1:]]
    ivar = [image[1].data["IVAR"][0, 1:]]
    
    p_opt, p_cov, continuum, rectified_model_flux  = model.fit(
        wl,
        flux,
        ivar,
        #R=4000,
        vacuum=False, 
        continuum_basis=Sinusoids(deg=3)
    )
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True)
    
    for i in range(len(wl)):
        axes[0].plot(wl[i], flux[i], c='k')
        axes[0].plot(wl[i], continuum[i], c="tab:red")
        axes[1].plot(wl[i], flux[i] / continuum[i], c='k')
        axes[1].plot(wl[i], rectified_model_flux[i], c="tab:red")    

    #fig.axes[0].set_title(f"index = {index}")   
    ''' 
    
    
    
    from specutils import Spectrum1D
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/j174239-133332red_multi.fits")# flux_ext=6)
    import matplotlib.pyplot as plt
    trim_blue, trim_red = (20, 20)
    #v_rel = -42.3
    v_rel = 35.0
    wl, flux, ivar = ([], [], [])
    for index in range(0, len(spectra)):
        wl.append(spectra[index].wavelength[trim_blue:-trim_red] * (1 - v_rel/299792.458))
        flux.append(spectra[index].flux[trim_blue:-trim_red])
        ivar.append(spectra[index].ivar[trim_blue:-trim_red])   

    blue_spectra = Spectrum1D.read("/Users/andycasey/Downloads/j174239-133332blue_multi.fits")# flux_ext=6)
    import matplotlib.pyplot as plt
    trim_blue, trim_red = (20, 20)
    for index in range(0, len(blue_spectra)):
        wl.append(blue_spectra[index].wavelength[trim_blue:-trim_red] * (1 - v_rel/299792.458))
        flux.append(blue_spectra[index].flux[trim_blue:-trim_red])
        ivar.append(blue_spectra[index].ivar[trim_blue:-trim_red])       
    '''
    p_opt, p_cov, continuum, rectified_model_flux = model.fit(
        wl, 
        flux, 
        ivar, 
        R=40_000,
        vacuum=False, 
        continuum_basis=Polynomial(deg=0),
    )

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2)# gridspec_kwds=dict(height_ratios=[2, 1]))
    
    for i in range(len(continuum)):
        axes[0].plot(wl[i], flux[i], c='k')
        axes[0].plot(wl[i], continuum[i], c="tab:red")
        axes[1].plot(wl[i], flux[i] / continuum[i], c='k')
        axes[1].plot(wl[i], rectified_model_flux[i], c="tab:red")
    
    '''
    deg = 3
    C, _ = model.basis_vectors.shape
    
    p_opt, p_cov, continuum, rectified_model_flux = model.fit(
        wl, 
        flux, 
        ivar, 
        R=40_000,
        vacuum=False, 
        continuum_basis=Sinusoids(deg=deg)
    )

    
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True)
    
    for i in range(len(wl)):
        axes[0].plot(wl[i], flux[i], c='k')
        axes[0].plot(wl[i], continuum[i], c="tab:red")
        axes[1].plot(wl[i], flux[i] / continuum[i], c='k')
        axes[1].plot(wl[i], rectified_model_flux[i], c="tab:red")    
    
    axes[1].set_ylim(0, 1.2)
    for ax in axes:
        ax.set_xlim(3000, 9000)
    
    raise a
        
    
    
    from specutils import Spectrum1D
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563blue_multi.fits")# flux_ext=6)
    import matplotlib.pyplot as plt
    trim_blue, trim_red = (20, 20)
    v_rel = -42.3
    wl, flux, ivar = ([], [], [])
    for index in range(0, len(spectra)):
        wl.append(spectra[index].wavelength[trim_blue:-trim_red] * (1 - v_rel/299792.458))
        flux.append(spectra[index].flux[trim_blue:-trim_red])
        ivar.append(spectra[index].ivar[trim_blue:-trim_red])   
    
    model = ContinuumModel.from_path("bosz-highres-optical-basis-vectors.h5")
    
    deg = 3
    C, _ = model.basis_vectors.shape
    
    p_opt, p_cov, continuum, rectified_model_flux = model.fit(
        wl, 
        flux, 
        ivar, 
        R=40_000,
        p0=p_opt[:32],
        vacuum=False, 
        continuum_basis=Sinusoids(deg=deg)
    )

    
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True)
    
    for i in range(len(wl)):
        axes[0].plot(wl[i], flux[i], c='k')
        axes[0].plot(wl[i], continuum[i], c="tab:red")
        axes[1].plot(wl[i], flux[i] / continuum[i], c='k')
        axes[1].plot(wl[i], rectified_model_flux[i], c="tab:red")    
    
    axes[1].set_ylim(0, 1.2)
    for ax in axes:
        ax.set_xlim(wl[0][0], wl[-1][-1])
            
            
    raise a
        
    for i in range(O[-1]):
        mask = (oi == i)
        #ls = ":" if fig is None else "-"
        ls = "-"
        mean = np.mean(wavelength[i])
        ax.plot(wls - mean, i + rectified_model_flux, c="tab:red", ls=ls)
        ax.plot(wavelength[i] - mean, i + flux[i] / continuum[mask], c='k', ls=ls)
        #sw, ew = wls[mask][[0, -1]]
        #ax.set_xlim(sw - 10, ew + 10)
        
    raise a
    
    '''
    A, basis_vectors, theta, flux, inv_sigma = model.fit(
        wl[:2], 
        flux[:2], 
        ivar[:2], 
        #R=40_000,
        vacuum=False, 
        continuum_basis=Sinusoids(deg=4),
    )
    
        
    C = 32
    
    def f(theta):
        return ((1 - theta[:C] @ basis_vectors) * (A @ theta[C:]) - flux) * inv_sigma
    
    def jacobian(theta):
        theta = np.atleast_2d(theta)
        return np.hstack([            
            (-basis_vectors * (theta[:, C:] @ A.T)).T,
            (1 - theta[:, :C] @ basis_vectors).T * A
        ]) * inv_sigma[:, None]
    
    theta = np.random.uniform(size=theta.shape)

    print("computing mine")
    mine = jacobian(theta)
    print("autograd")    

    jac = autograd_jacobian(f)
    truth = jac(theta)
    
    for index in range(theta.size):
        fig, ax = plt.subplots()
        ax.scatter(truth[:, index], mine[:, index])# - foo[:, index])
        ax.set_title(index)
    
    raise a
    '''
    
    
    
    '''
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563red_multi.fits", flux_ext=6)
    for index in range(0, len(spectra)):
        wl.append(spectra[index].wavelength[trim_blue:] * (1 - v_rel/299792.458))
        flux.append(spectra[index].flux[trim_blue:])
        ivar.append(spectra[index].ivar[trim_blue:])        
    '''

    raise a    
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