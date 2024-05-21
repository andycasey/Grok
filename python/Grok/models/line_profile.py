import numpy as np
from typing import Optional, Sequence, Union#, Self
from scipy import stats, optimize as op
from Grok.specutils.spectrum import Spectrum
from Grok.specutils.utils import ivar_to_sigma

from Grok.models.continuum_basis import ContinuumBasis


class BaseLineProfile:

    def __init__(
        self,
        λ: float,
        window: Optional[float] = 2,
        λ_tolerance: Optional[float] = 0.1,
        continuum_basis: Optional[ContinuumBasis] = None
    ):
        # TODO: add mask, other options, etc
        self.λ = λ
        self.window = window
        self.λ_tolerance = λ_tolerance
        self.continuum_basis = continuum_basis
        return None
        

    def fit(self, spectrum: Spectrum, initial: Optional[float] = None, op_kwds: Optional[dict] = None, **kwargs):
        
        # Get the data.
        si, ei = spectrum.λ.searchsorted([self.λ - self.window, self.λ + self.window]) + [0, 1]
        λ, flux, ivar = (spectrum.λ[si:ei], spectrum.flux[si:ei], spectrum.ivar[si:ei])
        sigma = ivar_to_sigma(ivar)
                    
        # Define a (faster) callable function to optimize (which should match __call__!), and set initial values
        f, p0 = self.get_callable_and_initial_guess(λ, flux, ivar)
                        
        kwds = dict(xtol=1e-12, ftol=1e-12, verbose=0)
        kwds.update(op_kwds or {})                        
            
        popt, pcov, infodict, mesg, ier = op.curve_fit(
            f,
            λ,
            flux,
            p0=p0,
            sigma=sigma,
            bounds=self.get_bounds(**kwargs),
            full_output=True,
            **kwds
        )
                
        self.theta_, self.theta_cov_ = (popt, pcov)
        # TODO: do checks on quality of the fit, etc.. should we raise an exception if the fit is bad?        
        return self


    def __call__(self, λ, *θ):                    
        if self.continuum_basis is None:
            return self._profile(λ, *θ)
        else:
            A = self.continuum_basis.design_matrix(λ)
            return (A @ θ[self.num_profile_parameters:]) * self._profile(λ, *θ)
        
    
        
class GaussianLineProfile(BaseLineProfile):
    
    num_profile_parameters = 3
    
    def get_callable_and_initial_guess(self, λ, flux, ivar):        
        p0 = [self.λ, 0.1, 0.2]
        if self.continuum_basis is not None:
            A = self.continuum_basis.design_matrix(λ)
            # TODO: consider different logic for initial guess if we have a NMF model                        
            p0.extend(self.continuum_basis.get_initial_guess(A, flux, ivar))
            def f(λ, *θ):
                return (A @ θ[self.num_profile_parameters:]) * self._profile(λ, *θ)            
        else:
            f = self._profile                
        return (f, p0)
            
    
    def get_bounds(self, σ_max=10, amplitude_max=1, **kwargs):
        bounds = [
            (self.λ - self.λ_tolerance, self.λ + self.λ_tolerance),
            (0, σ_max),
            (0, amplitude_max)
        ]
        if self.continuum_basis is not None:
            bounds.extend([(-np.inf, +np.inf)] * self.continuum_basis.num_parameters)
        return np.atleast_2d(bounds).T
        
    def _profile(self, λ, μ, σ, amplitude, *_, **__):
        return 1 - amplitude * stats.norm.pdf(λ, loc=μ, scale=σ)
        
    @property
    def ew(self):
        """Return the equivalent width of the fitted profile in milliAngstroms (mÅ)."""
        try:
            # μ, σ, amplitude, *p = self.theta_
            return 1000 * self.theta_[2]             
        except AttributeError:
            return None

    @property
    def u_ew(self):
        """Return the uncertainty on the equivalent width of the fitted profile in milliAngstroms (mÅ)."""
        return 1000 * self.theta_cov_[2, 2]**0.5        
    


# THINKOS:
# - VoigtLineProfile sub-class
# - can we avoid having `num_profile_parameters` as a class attribute?



if __name__ == "__main__":
    
    
    from Grok.specutils import Spectrum
    from Grok.models.continuum_basis import PolynomialBasis, FourierBasis
    
    spectrum = Spectrum.read("/Users/andycasey/software/smhr/smh/data/spectra/hd122563.fits")
    
    f = GaussianLineProfile(λ=4736.773, continuum_basis=PolynomialBasis(1)).fit(spectrum)
    
    import matplotlib.pyplot as plt
    fig = spectrum.plot()
    ax = fig.axes[0]
    ax.plot(spectrum.λ, f(spectrum.λ, *f.theta_), c='r')
    ax.set_xlim(4736.773 - f.window, 4736.773 + f.window)
    ax.set_ylim(0, 1.2)
    