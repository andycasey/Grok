
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
from astropy.constants import c

from Grok.specutils.lsf import instrument_lsf_kernel
from Grok.specutils.spectrum import Spectrum, apply_relativistic_doppler_shift

from threading import Thread

from typing import Iterable

C_KM_S = c.to("km/s").value

def apply_radial_velocity_shift(λ, v_rel):
    beta = v_rel / C_KM_S
    scale = np.sqrt((1 + beta) / (1 - beta))     
    return λ * scale

class cumulative:
    
    """A cumulative iterator."""
    
    def __init__(self, iterable: Iterable, start: Optional[int] = 0):
        self.iterable = iter(iterable)
        self.sum = start
        return None
    
    def __iter__(self):
        return self
    
    def __next__(self):
        item = next(self.iterable)
        I = len(item)
        self.sum += I
        return (self.sum - I, item)        


def expand_path(path):
    return os.path.abspath(os.path.expanduser(path))




def get_L_default(wavelength):
    return 2 * np.ptp(wavelength)


class ContinuumBasis:
    pass

class NoContinuum(ContinuumBasis):
    
    @property
    def num_parameters(self):
        return 0

    def design_matrix(self, wavelength):
        return 0

class Sinusoids(ContinuumBasis):    
        
    def __init__(self, P=7, L=None):
        self.P = P
        self.L = L
        return None
    
    @property
    def num_parameters(self):
        return self.P
        
    def default_L(self, wavelength):
        return 2 * np.ptp(wavelength)
    
    def design_matrix(self, wavelength):
        if self.L is None:
            L = self.default_L(wavelength)
        elif isinstance(self.L, (float, int)):
            L = self.L
        else:
            L = self.L(wavelength)
            
        scale = (np.pi * wavelength) / L
        A = np.ones((wavelength.size, self.P), dtype=float)
        for j, f in zip(range(1, self.P), cycle((np.sin, np.cos))):
            A[:, j] = f(scale * (j + (j % 2)))        
        return A

            
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
    
    def __init__(
        self, 
        basis_vacuum_wavelength, 
        stellar_basis_vectors, 
        telluric_basis_vectors=None,
        meta=None,         
    ):
        """
        :param basis_vacuum_wavelength:
            An array of vacuum wavelengths (in Angstroms) of the basis vectors.
        
        :param stellar_basis_vectors:
            A (C, P)-shape array of stellar basis vectors with non-negative entries. Here, C is
            the number of basis vectors and P is the same as the size of `basis_vacuum_wavelength`.
                    
        :param telluric_basis_vectors: [optional]
            A (B, P)-shape array of telluric basis vectors with non-negative entries. Here, P is
            the same as the size of `basis_vacuum_wavelength`.
            
        :param meta: [optional]
            A metadata dictionary.            
        """
        self.basis_vacuum_wavelength = basis_vacuum_wavelength
        self.stellar_basis_vectors = stellar_basis_vectors            
        self.telluric_basis_vectors = telluric_basis_vectors
        self.meta = meta or dict()   
        return None        
    

    @property
    def can_model_tellurics(self):
        return self.telluric_basis_vectors is not None 
            
            
    @classmethod
    def from_path(cls, path):
        """
        Load a continuum model from a file.
        
        :param path:
            The path to the continuum model file.
        """
        with h5.File(expand_path(path), "r") as fp:
            vacuum_wavelength = fp["vacuum_wavelength"][:]
            stellar_basis_vectors = fp["stellar_basis_vectors"][:]
            meta = dict(fp["stellar_basis_vectors"].attrs)       
            try:
                telluric_basis_vectors = fp["telluric_basis_vectors"][:]
            except:
                telluric_basis_vectors = None
                
        return cls(
            vacuum_wavelength,
            stellar_basis_vectors,
            telluric_basis_vectors=telluric_basis_vectors,
            meta=meta,            
        )
    
    
    def resample_basis_vectors(self, λ_vacuum, stellar_v_rel=0.0, telluric_v_rel=0.0, Ro=None, Ri=None, telluric=True, **kwargs):
        """
        Resample (and optionally convolve) the basis vectors to the observed wavelengths.
        
        :param λ_vacuum:
            The observed vacuum wavelengths.
        
        :param Ro: [optional]
            The spectral resolution of the observations. If `None` is given, then no
            convolution will take place; the basis vectors will be interpolated.
        
        :param Ri: [optional]
            The spectral resolution of the basis vectors. If `None` is given then
            it defaults to the input spectral resolution stored in the metadata of
            the basis vector file, or infinity if the input spectral resolution is
            not stored in that file.
                
        """

        λ_vacuum_stellar = apply_radial_velocity_shift(λ_vacuum, stellar_v_rel)
        λ_vacuum_telluric = apply_relativistic_doppler_shift(λ_vacuum, telluric_v_rel)
        
        
        if Ro is None:
            # Interpolation only.            
                
            basis_vectors = [_interpolate_basis_vector(λ_vacuum_stellar, self.basis_vacuum_wavelength, self.stellar_basis_vectors)]
            if telluric:
                basis_vectors.append(_interpolate_basis_vector(λ_vacuum_telluric, self.basis_vacuum_wavelength, self.telluric_basis_vectors))
            
        else:
            # Convolution
            Ri = Ri or self.meta.get("Ri", np.inf)
            
            # special case to avoid double-building the same convolution kernel            
            K = instrument_lsf_kernel(self.basis_vacuum_wavelength, λ_vacuum_stellar, Ro, Ri)                
            
            basis_vectors = [self.stellar_basis_vectors @ K]
            if telluric:
                # Only reconstruct the kernel if the v_rel != telluric_v_rel
                if stellar_v_rel != telluric_v_rel:
                    K = instrument_lsf_kernel(self.basis_vacuum_wavelength, λ_vacuum_telluric, Ro, Ri)
                
                basis_vectors.append(self.telluric_basis_vectors @ K)
    
        basis_vectors = np.vstack(basis_vectors)
        S = self.stellar_basis_vectors.shape[0]
        T = basis_vectors.shape[0] - S

        return (basis_vectors, S, T)
    


    def get_initial_guess(self, A, flux, ivar, order_masks, n_continuum_params, basis_vectors, p0=None, small=1e-6):
        n_bv = basis_vectors.shape[0]
        x0 = small * np.ones(n_bv + A.shape[1])
        if p0 is None:
            scalar = np.ones_like(flux)
        else:
            scalar = 1 - p0 @ basis_vectors
            x0[:n_bv] = p0
        si = 0
        for n, mask in zip(n_continuum_params, order_masks):      
            scale = scalar[mask]      
            x0[n_bv + si:n_bv + si + n] = _solve_X(flux[mask] / scale, ivar[mask] * scale**2, A[mask, si:si+n])
            si += n
        return x0
    


    
    def fit(
        self,
        spectra: Sequence[Spectrum],
        continuum_basis: Union[ContinuumBasis, Sequence[ContinuumBasis]] = Sinusoids,
        p0: Optional[Sequence[float]] = None,
        stellar_v_rel: Optional[float] = 0.0,
        telluric_v_rel: Optional[float] = 0.0,
        R: Optional[float] = None,
        telluric: Optional[bool] = True,
        callback: Optional[callable] = None,
        op_kwds: Optional[dict] = dict(xtol=1e-16, verbose=2),
        **kwargs
    ):
        """
        Fit the model to some spectra in the observed frame.
        
        :param spectra:
            A sequence of observed spectra.
        
        :param continuum_basis: [optional]
            The continuum basis to use. This can be a single instance of a continuum basis
            or a sequence of continuum bases. If a sequence is given, then the length must
            match the number of orders in the spectra.
        
        :param R: [optional]
            The spectral resolution of the observations. If `None` is given, then no
            convolution will take place; the basis vectors will be interpolated.
        
        :param p0: [optional]
            An initial guess for the model parameters.
        
        :param telluric: [optional]
            Specify whether to model tellurics.
        
        :param callback: [optional]
            A function to call after each iteration of the optimization.
        """
        
        λ_vacuum, flux, ivar, oi, pi, S, P = _prepare_spectra(spectra)
                
        # Expand the continuum bases.
        continuum_basis = _expand_continuum_basis(continuum_basis, S)
                    
        # Resample the stellar and telluric basis vectors, if the v_rel has changed.
        
        #try:
        #    basis_vectors, n_stellar_bases, n_telluric_bases = self.resampled_basis_vectors
        #except AttributeError:            
        basis_vectors, n_stellar_bases, n_telluric_bases = self.resampled_basis_vectors = self.resample_basis_vectors(
            λ_vacuum, 
            stellar_v_rel=stellar_v_rel, 
            telluric_v_rel=telluric_v_rel,
            Ro=R, 
            telluric=telluric,
            **kwargs
        )
        
        # Tellurics are usually less flexibly modelled (there are fewer basis vectors than
        # the stellar basis vectors), so we need to mask them out.
        if telluric:
            telluric_mask_min = kwargs.pop("telluric_mask_min", 0.01)
            telluric_mask = np.any(basis_vectors[n_stellar_bases:] >= telluric_mask_min, axis=0)
            ivar[telluric_mask] *= 1e-4
        
        # Construct the packed design matrix.
        n_continuum_params = [b.num_parameters for b in continuum_basis]
                        
        A = np.zeros((λ_vacuum.size, sum(n_continuum_params)))
        si, order_masks = (0, [])
        for o, (cb, n) in enumerate(zip(continuum_basis, n_continuum_params)):
            mask = (oi == o)
            order_masks.append(mask)
            A[mask, si:si+n] = cb.design_matrix(λ_vacuum[mask])
            si += n

        # Get an initial guess of continuum with no absorption.
        x0 = self.get_initial_guess(A, flux, ivar, order_masks, n_continuum_params, basis_vectors, p0)
        
        # Let's pre-compute some things for faster Jacobian evaluations.
        n_bases = n_stellar_bases + n_telluric_bases
        inv_sigma = np.sqrt(ivar)        
        mBVT = -basis_vectors.T
        A_inv_sigma = A * inv_sigma[:, None]
        mBVT_inv_sigma = mBVT * inv_sigma[:, None]
        
        if sum(n_continuum_params) > 0:            
            def jacobian(theta):
                return np.hstack([
                    (A @ theta[n_bases:, None]) * mBVT_inv_sigma,
                    (1 + mBVT @ theta[:n_bases, None].T) * A_inv_sigma
                ])

            def f(theta):
                y = (1 - theta[:n_bases] @ basis_vectors) 
                if A.shape[1] > 0:
                    y *= (A @ theta[n_bases:])
                return (y - flux) * inv_sigma         

        else:
            def f(theta):
                return (1 - theta @ basis_vectors - flux) * inv_sigma                        

            def jacobian(x):
                return mBVT_inv_sigma
                        
        print("Starting")
        t_init = time()
        result = op.least_squares(
            f, 
            x0, 
            jac=jacobian, 
            bounds=self.get_bounds(x0.size, n_bases),
            **op_kwds
        )
        t_op = time() - t_init
        print(f"took {t_op:.1f} s to optimize {len(x0)} parameters with {flux.size} data points")
                                
        if sum(n_continuum_params) > 0:       
            y_pred = [A @ result.x[n_bases:]] # continuum
        else:
            y_pred = [np.ones_like(flux)]
        
        y_pred.append(       
            1 - result.x[:n_stellar_bases] @ basis_vectors[:n_stellar_bases] # rectified_stellar_flux
        )
        if n_telluric_bases > 0:
            # rectified_telluric_flux
            y_pred.append(
                1 - result.x[n_stellar_bases:n_bases] @ basis_vectors[n_stellar_bases:n_bases],
            )
        else:
            y_pred.append(np.nan * np.ones_like(y_pred[0])) # rectified_telluric_flux        
                
        y_pred = np.array(y_pred)
        
        continuum, rectified_model_flux, rectified_telluric_flux = zip(*[y_pred[:, m] for m in order_masks])
        return (result, continuum, rectified_model_flux, rectified_telluric_flux)
                    

    def get_bounds(self, n_params, n_basis_vectors, component_bounds=(1e-10, +np.inf)):
        return np.vstack([
            np.tile(component_bounds, n_basis_vectors).reshape((n_basis_vectors, 2)),
            np.tile([-np.inf, +np.inf], n_params - n_basis_vectors).reshape((-1, 2))
        ]).T            


def _solve_X(Y: Sequence[float], Cinv: Sequence[float], A: np.array):
    MTM = A.T @ (Cinv[:, None] * A)
    MTy = A.T @ (Cinv * Y)
    return np.linalg.solve(MTM, MTy)
        
        
def instantiate(item, **kwargs):
    if isinstance(item, type):
        return item(**kwargs)
    else:
        return item

def _prepare_spectra(spectra):
    spectra = [spectra] if isinstance(spectra, Spectrum) else spectra
    S = len(spectra)
    P = sum(tuple(map(len, spectra)))
    
    λ_vacuum, flux, ivar, oi = (np.empty(P), np.empty(P), np.empty(P), np.empty(P, dtype=int))
    for i, (si, spectrum) in enumerate(cumulative(spectra)):
        ei = si + spectrum.λ.size
        oi[si:ei] = i
        λ_vacuum[si:ei] = spectrum.λ_vacuum
        flux[si:ei] = spectrum.flux
        ivar[si:ei] = spectrum.ivar

    pi = np.argsort(λ_vacuum) # pixel indices
    λ_vacuum, flux, ivar, oi = (λ_vacuum[pi], flux[pi], ivar[pi], oi[pi])
    
    bad_pixel = (~np.isfinite(flux)) | (~np.isfinite(ivar))
    ivar[bad_pixel] = 0
    return (λ_vacuum, flux, ivar, oi, pi, S, P)


def _expand_continuum_basis(continuum_basis, S):
    if isinstance(continuum_basis, (list, tuple)):
        if len(continuum_basis) != S:
            raise ValueError(f"a sequence of continuum_basis was given, but the length does not match the number of spectra ({len(continuum_basis)} != {S})")
        return tuple(map(instantiate, continuum_basis))
    else:
        return tuple([instantiate(continuum_basis)] * S)




def _interpolate_basis_vector(λ, basis_vacuum_wavelength, basis_vectors):
    C = basis_vectors.shape[0]
    P = λ.size
    bv = np.zeros((C, P))
    for c, basis_vector in enumerate(basis_vectors):
        bv[c] = np.interp(λ, basis_vacuum_wavelength, basis_vector, left=0, right=0)
    return bv



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
    
    try:
        model
    except NameError:
        model = ContinuumModel.from_path("Grok/bosz-highres-optical-basis-vectors.h5")
        #from Grok.specutils_new import Spectrum
        '''
        orig_spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563red_multi.fits")# flux_ext=6)
        import matplotlib.pyplot as plt
        trim_blue, trim_red = (20, 20)
        spectra = []
        for index in range(0, len(orig_spectra) - 2):
            wl = orig_spectra[index].wavelength[trim_blue:-trim_red] #* (1 - v_rel/299792.458)
            flux = orig_spectra[index].flux[trim_blue:-trim_red]
            ivar = orig_spectra[index].ivar[trim_blue:-trim_red]
            
            spectra.append(Spectrum(wl, flux, ivar, medium="air"))
        '''
        '''
        spectra = []
        from astropy.io import fits
        with fits.open("/Users/andycasey/Downloads/pepsir.20230619.018.dxt.nor") as image:
            wl = image[1].data["Arg"]
            flux = image[1].data["Fun"]
            ivar = 1.0/image[1].data["Var"]
            spectra.append(Spectrum(wl, flux, ivar, medium="air"))
        '''
        
    else:
        print("Using previous model and spectra!!!")
    
    from Grok.specutils.spectrum import Spectrum, SpectrumCollection
    
    spectra = [
        Spectrum.read("/Users/andycasey/Downloads/pepsir.20230619.018.dxt.nor")
    ]
    spectra = SpectrumCollection.read("/Users/andycasey/Downloads/hd122563red_multi.fits")
    
    #v_rel = -20166.342/1000.0
    #spectra[0].apply_velocity_shift(v_rel)
    
    
    (result, continuum, rectified_model_flux, rectified_telluric_flux) = model.fit(
        spectra, 
        stellar_v_rel=0.0,
        #telluric_v_rel=-20166.342/1000.0,
        telluric=True,
        continuum_basis=NoContinuum(),
        xtol=1e-16,
        ftol=1e-64,
    )
    
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True)
    
    for i, spectrum in enumerate(spectra):
        
        axes[0].plot(spectrum.λ, spectrum.flux, c='k', label="Data" if i == 0 else None)
        axes[0].plot(spectrum.λ, continuum[i] * rectified_model_flux[i], c="tab:red", label="Stellar absorption" if i == 0 else None)    
        axes[0].plot(spectrum.λ, continuum[i] * rectified_telluric_flux[i], c="tab:blue", label="Telluric" if i == 0 else None)    
        #axes[0].plot(spectrum.λ, continuum[i], c="#666666")
        
        axes[1].plot(spectrum.λ, spectrum.flux / continuum[i], c='k')
        axes[1].plot(spectrum.λ, rectified_model_flux[i], c="tab:red")    
        axes[1].plot(spectrum.λ, rectified_telluric_flux[i], c="tab:blue")
    
    axes[0].legend(frameon=False)
    axes[1].set_ylim(0, 1.2)
    axes[0].set_ylim(0, 1.2)
    axes[1].set_xlabel("Wavelength [A]")
    axes[0].set_xlim(spectra[0].λ[0], spectra[-1].λ[-1])
    
    raise a

    # Given some line list, we want to:
    # - fit some profile at each wavelength;
    # - use the existing stellar model as the initial guess of the continuum;
    # - use the existing stellar model as the initial guess of the line profile;
    # - generate masks in the windowed region of the line based on some heuristic.
    
    line_wavelengths = np.loadtxt("/Users/andycasey/Downloads/pepsi_linelist.moog", usecols=(0, ))
    
    window = 2 # angstroms
    
    min_rectified_model_stellar_flux = 0.995
    min_rectified_model_telluric_flux = 0.995
    
    spectrum, = spectra
    rectified_telluric_flux, = rectified_telluric_flux
    rectified_model_flux, = rectified_model_flux
    
    continuum_basis = None
    
    def to_contiguous_regions(mask):
        v = np.diff(mask.astype(int))
        indices = np.hstack([
            0, 
            np.repeat(1 + np.where(v != 0)[0], 2),
            len(mask)
        ])
                
        indices = indices.reshape((-1, 2))
        if indices.size > 0:
            offset = 0 if mask[0] else 1
            return indices[offset::2]
        return indices
            
    
    

    # For each line, we want to fit a Gaussian profile.
    for λ in air_to_vacuum(line_wavelengths):        
        if not (spectrum.λ[0] < λ < spectrum.λ[-1]):
            continue
            
        # use the existing spectral model.
        si, ei = spectrum.λ.searchsorted([λ - window, λ + window]) + [0, 1]
        
        Y = spectrum.flux[si:ei]
        Cinv = spectrum.ivar[si:ei]
        
        # generate a mask
        mask = (
            (rectified_telluric_flux[si:ei] < min_rectified_model_telluric_flux)
        |   (rectified_model_flux[si:ei] < min_rectified_model_stellar_flux)
        )
        
        fig, ax = plt.subplots()
        ax.plot(spectrum.λ[si:ei], Y, c='k')
        ax.plot(spectrum.λ[si:ei], rectified_model_flux[si:ei], c="tab:red")
        ax.plot(spectrum.λ[si:ei], rectified_telluric_flux[si:ei], c="tab:blue")
        
        for sm, em in to_contiguous_regions(mask):
            ax.axvspan(spectrum.λ[si + sm], spectrum.λ[si + em], color="#cccccc", alpha=0.2)
        
        ax.axvline(λ, c="#666666", ls=":")
        
        
    
    
    
    
    from specutils import Spectrum1D
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563red_multi.fits")# flux_ext=6)
    import matplotlib.pyplot as plt
    trim_blue, trim_red = (20, 20)
    v_rel = -42.3
    #v_rel = 35.0
    v_rel = 69.1
    wl, flux, ivar = ([], [], [])
    for index in range(0, len(spectra)):
        wl.append(spectra[index].wavelength[trim_blue:-trim_red] * (1 - v_rel/299792.458))
        flux.append(spectra[index].flux[trim_blue:-trim_red])
        ivar.append(spectra[index].ivar[trim_blue:-trim_red])   

    blue_spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563blue_multi.fits")# flux_ext=6)
    import matplotlib.pyplot as plt
    trim_blue, trim_red = (20, 20)
    for index in range(0, len(blue_spectra)):
        wl.append(blue_spectra[index].wavelength[trim_blue:-trim_red] * (1 - v_rel/299792.458))
        flux.append(blue_spectra[index].flux[trim_blue:-trim_red])
        ivar.append(blue_spectra[index].ivar[trim_blue:-trim_red])
    
    '''
    result, continuum, rectified_model_flux = model.fit(
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
    
    result, continuum, rectified_model_flux = model.fit(
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

    telluric = np.loadtxt("/Users/andycasey/Downloads/input/spectra/templates/Synth.Tellurics.350_1100nm/template.txt", skiprows=1)    
    axes[1].plot(telluric.T[0] * 10 * (1 - v_rel/299792.458), telluric.T[1], c="tab:blue")
    
    raise a
        
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True)
    
    axes[0].plot(wl, fl, c='k')
    axes[0].plot(wl, continuum[0], c="tab:red")
    axes[1].plot(wl, fl / continuum[0], c='k')
    axes[1].plot(wl, rectified_model_flux[0], c="tab:red")           

    
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
    
    result, continuum, rectified_model_flux = model.fit(
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