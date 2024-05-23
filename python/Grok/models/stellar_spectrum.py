
import numpy as np
from scipy import optimize as op
from types import MappingProxyType
from typing import Sequence, Union, Optional, Tuple

from Grok.models.basis import LinearBasis
from Grok.specutils.utils import ivar_to_sigma
from Grok.utils import overlap

class LinearStellarSpectrumModel:
    
    """Model a stellar spectrum as a (nearly) linear combination of basis spectra."""
    
    def __init__(
        self,
        spectral_basis: LinearBasis,
        continuum_basis: Optional[LinearBasis] = None,
        mask_regions: Optional[Sequence[Tuple[float, float]]] = None,
    ):
        """
        Model a stellar spectrum as a (nearly) linear combination of basis spectra.
        
        :param spectral_basis:
            The basis of spectral components. For example, a `Grok.models.NMFSpectralBasis` instance.
            
        :param continuum_basis: [optional]
            The basis of continuum components.
        """
        self.bases = MappingProxyType(dict(
            [(name, base) for name, base in zip(("spectral", "continuum"), (spectral_basis, continuum_basis)) if base is not None]            
        ))
        self.mask_regions = mask_regions
        return None
    
    
    @property
    def n_parameters(self):
        """The total number of model parameters."""
        return sum(tuple(base.n_parameters for base in self.bases.values()))
    
    def get_mask(self, λ):
        """
        Get a mask for the regions to exclude from the fit.
        
        :param λ:
            The wavelengths to fit.
        """
        mask = np.zeros(λ.size, dtype=bool)
        if self.mask_regions is not None:
            for region in self.mask_regions:
                mask |= (region[1] >= λ) & (λ >= region[0])
        return mask

    
    def prepare_fit(self, spectrum, **kwargs):
        """
        Given a spectrum and a window to fit around, prepare the data and design matrices.
        
        :param spectrum:
            The spectrum to fit.
            
        :param window:
            The window to fit around the central wavelength. The full fitting domain is twice this window.
        """
        
        print("Change this to allow for `spectra`, not just `spectrum`.")
        
        # Get the data.
        λ_min, λ_max = overlap(self.bases["spectral"].λ_rest_vacuum, spectrum.λ)
        
        si, ei = spectrum.λ.searchsorted([λ_min, λ_max]) + [0, 1]        
        λ, flux, ivar = (spectrum.λ[si:ei], spectrum.flux[si:ei], spectrum.ivar[si:ei])
    
        use = ~self.get_mask(λ)
            
        # Build design matrices and get bounds.

        bounds = self.bases["spectral"].get_bounds()
        # TODO: this could cause unnecessary edge effects by masking the λ beforehand
        As = self.bases["spectral"].get_design_matrix(λ[use], **kwargs)
        
        if "continuum" in self.bases:
            Ac = self.bases["continuum"].get_design_matrix(λ[use])
            bounds.extend(self.bases["continuum"].get_bounds())
        else:
            Ac = None
                
        return (λ, flux, ivar, use, As, Ac, bounds)        


            
    def fit(
        self,
        spectrum,
        p0: Optional[Sequence[float]] = None,
        stellar_v_rel: Optional[float] = None,
        telluric_v_rel: Optional[float] = None,
        Ro: Optional[float] = None,
        **kwargs
    ):
        """
        Fit a linear stellar spectrum model to the given data.
        
        :param spectrum:
            The spectrum to fit.
        
        :param stellar_v_rel: [optional]
            The relative velocity of the star.
        
        :param telluric_v_rel: [optional]
            The relative velocity of the telluric lines.
        
        :param Ro: [optional]
            The resolution of the observed spectrum.        
        """
        
        (λ, flux, ivar, use, As, Ac, bounds) = self.prepare_fit(
            spectrum,
            stellar_v_rel=stellar_v_rel,
            telluric_v_rel=telluric_v_rel,
            Ro=Ro,
            **kwargs
        )
            
        n = self.bases["spectral"].n_parameters
            
        if Ac is None:
            def f(λ, *θ):
                return 1 - As @ θ # TODO: this assumes the stellar basis is an absorption type (e.g., NMF)
        else:
            def f(λ, *θ):
                return (Ac @ θ[n:]) * (1 - As @ θ[:n])

        if p0 is None:
            #p0 = np.mean(bounds, axis=1)
            #p0[~np.isfinite(p0)] = 1
            p0 = 1e-5 * np.ones(self.n_parameters)
            p0[n:] = 0

        self.θ, self.Σ = op.curve_fit(
            f,
            λ[use],
            flux[use],
            p0=p0,
            sigma=ivar_to_sigma(ivar[use]),
            bounds=np.array(bounds).T
        )
        
        return (λ[use], f(λ[use], *self.θ))
    
    
    
if __name__ == "__main__":


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

    
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    
    from Grok.models.basis import NMFSpectralBasis, PolynomialBasis, FourierBasis
    
    spectral_basis = NMFSpectralBasis.from_path("Grok/bosz-highres-optical-basis-vectors.h5")    
    continuum_basis = PolynomialBasis(1)
    continuum_basis = FourierBasis(5)
    
    edge = 50
    c = 5891.5 #  air_to_vacuum(5655.493)
    mask_regions = [
        (0, c - edge),
        (c + edge, 1e4)
    ]
    
    model = LinearStellarSpectrumModel(spectral_basis, continuum_basis, mask_regions)
        
    from Grok.specutils import Spectrum
    
    spectrum = Spectrum.read("/Users/andycasey/software/smhr/smh/data/spectra/hd122563.fits")
    spectrum.λ = air_to_vacuum(spectrum.λ)
    
    λ_pred, y_pred = model.fit(spectrum, Ro=27_000)
    
    fig, ax = plt.subplots()
    ax.plot(spectrum.λ, spectrum.flux, c="k", label="Data")
    ax.plot(λ_pred, y_pred, c="r", label="Model")
    ax.set_xlim(λ_pred[[0, -1]])
    ax.set_ylim(0, 1.2)
    ax.axvline(c, c="#666666", ls=":")