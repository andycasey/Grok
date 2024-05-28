
import numpy as np
from scipy import optimize as op
from types import MappingProxyType
from typing import Sequence, Dict, Union, Optional, Tuple

from Grok.models.basis import LinearBasis
from Grok.specutils.utils import ivar_to_sigma
from Grok.utils import overlap

class LinearStellarSpectrumModel:
    
    """Model a stellar spectrum as a (nearly) linear combination of basis spectra."""
    
    def __init__(
        self,
        bases: Dict[str, LinearBasis],
        mask_regions: Optional[Sequence[Tuple[float, float]]] = None,
    ):
        """
        Model a stellar spectrum as a (nearly) linear combination of basis spectra.
        
        :param bases:
            A dictionary of linear bases. The keys are the names of the dictionary. This can include
            things like absorption spectra, telluric features, and continuum bases.
        """
        self.bases = MappingProxyType(dict([(k, v) for k, v in bases.items() if v is not None]))
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
    
    def get_bounds(self):
        return np.vstack([base.get_bounds() for base in self.bases.values()])
    
    def get_initial_guess(self, *args, **kwargs):
        p0 = []
        for base in self.bases.values():
            p0.extend(base.get_initial_guess(*args, **kwargs))
        return p0
    
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
        λ_min, λ_max = overlap(self.bases["stellar"].λ_rest_vacuum, spectrum.λ_rest_vacuum)
        si, ei = spectrum.λ_rest_vacuum.searchsorted([λ_min, λ_max]) + [0, 1]        
        λ, flux, ivar = (spectrum.λ_rest_vacuum[si:ei], spectrum.flux[si:ei], spectrum.ivar[si:ei])
    
        use = ~self.get_mask(λ)
            
        # Build design matrices and get bounds.
        p0, bounds, As = ([], [], [])
        for name, base in self.bases.items():
            A = base.get_design_matrix(λ[use], **kwargs)
            As.append(A)
            p0.extend(base.get_initial_guess(λ[use], flux[use], ivar[use]))
            bounds.extend(base.get_bounds())
                        
        return (λ, flux, ivar, use, As, bounds, p0)        

    
    def fit(self, spectrum, p0: Optional[Sequence[float]] = None, **kwargs):
        """
        Fit a linear stellar spectrum model to the given data.
        
        :param spectrum:
            The spectrum to fit.
        
        """
        
        (λ, flux, ivar, use, As, bounds, x0) = self.prepare_fit(spectrum, **kwargs)

        if p0 is None:
            p0 = x0
                    
        def f(λ, *θ, full_output=False):
            si, y = (0, 1)
            y_bases = {}
            for A, (name, base) in zip(As, self.bases.items()):
                y_bases[name] = base(λ, θ[si:si+A.shape[1]], A=A)
                y *= y_bases[name]
                si += A.shape[1]
            
            if full_output:
                return (y, y_bases)
            return y

        self.θ, self.Σ = op.curve_fit(
            f,
            λ[use],
            flux[use],
            p0=p0,
            sigma=ivar_to_sigma(ivar[use]),
            bounds=np.array(bounds).T
        )
        
        # Assign fitted variables to the bases
        si = 0
        for A, base in zip(As, self.bases.values()):
            n = A.shape[1]
            base.θ = self.θ[si:si+n]
            base.Σ = self.Σ[si:si+n, si:si+n]
            si += n
        
        print("save most recent design matrices")
        return self 
    

    
    def __call__(self, λ, *θ, **kwargs):
        
        si, y = (0, 1)
        y_bases = {}
        for n, A, (name, base) in zip(n_parameters, As, self.bases.items()):
            y_bases[name] = base(λ, θ[si:si+n], A=A)
            y *= y_bases[name]                
            si += n
        
        return y
        
    
        
    
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

        
    from Grok.specutils import Spectrum        
    from Grok.models.basis import NMFSpectralBasis, PolynomialBasis, FourierBasis
    from Grok.models.profile import VoigtProfile
    from Grok.models.spectral_line import SpectralLineModel
    from Grok.utils import overlap, to_contiguous_regions
    from copy import copy

    
    import h5py as h5
    with h5.File("Grok/bosz-highres-optical-basis-vectors.h5") as fp:
        stellar_base = NMFSpectralBasis(fp["vacuum_wavelength"][:], fp["stellar_basis_vectors"][:].T)
        telluric_base = NMFSpectralBasis(fp["vacuum_wavelength"][:], fp["telluric_basis_vectors"][:].T, v_rel=-5)
        

    ll = air_to_vacuum(np.loadtxt("/Users/andycasey/Downloads/pepsi_linelist.moog", usecols=(0, )))
    ll = air_to_vacuum(np.loadtxt("/Users/andycasey/Downloads/linelist_mm.txt", usecols=(0, )))

    #spectrum = Spectrum.read("/Users/andycasey/Downloads/pepsir.20230619.018.dxt.nor")
    #spectrum = Spectrum.read("/Users/andycasey/Downloads/pepsir.20231102.014.sxt.avr")
    spectrum = Spectrum.read("/Users/andycasey//Downloads/pepsib.20161117.040.dxt.all6.gz") # Solar spectrum?
    spectrum.apply_velocity_shift(-5)
    
    lower, upper = overlap(spectrum.λ_vacuum, ll)    
    #print(lower, upper)

    # Mask everything outside the overlap region.
    mask_regions = [(0, lower - 3), (upper + 3, 1e4)]
    
    bases = dict(
        stellar=stellar_base,
        telluric=telluric_base,
        #continuum=PolynomialBasis(0)
    )
    
    model = LinearStellarSpectrumModel(bases, mask_regions=mask_regions)
    
    λ_pred, (y_pred, y_bases) = model.fit(spectrum)

    # Note: make sure we are plotting the spectrum in vacuum frame, because that what the LinearStellarSpectrumModel fits in
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()        
    ax.plot(spectrum.λ_vacuum, spectrum.flux, c="k", label="Data")
    
    for name, base in model.bases.items():
        ax.plot(λ_pred, y_bases[name], label=name)
        
    ax.plot(λ_pred, y_pred, c="r", label="Model")    
    ax.set_xlim(λ_pred[[0, -1]])
    ax.set_ylim(0, 1.2)
    ax.legend()
    

    print(lower, upper)
    
    # Fit some EWs
    count = 0
    from time import time
    for λ in np.sort(ll):
        
        if not (upper > λ > lower):
            continue
            

        l, u = (λ - 0.15, λ + 0.15)
        bases = dict(
            stellar=copy(model.bases["stellar"]).set_taper_region(l, u, 0.025),
            telluric=copy(model.bases["telluric"]).set_taper_region(l, u, 0.025),
        )
        
        window = 1
        
        line_model = SpectralLineModel(λ, bases=bases)
        
        line_model.fit(spectrum, window=window)
                
        fig = line_model.plot(spectrum, window=window)
        
        '''
        lm2 = SpectralLineModel(λ, bases=bases)
        x0 = np.mean(lm2.get_bounds(), axis=1)
        x0[~np.isfinite(x0)] = 1e-5        
        lm2.fit(spectrum, window=window, p0=x0)
        
        lm2.plot(spectrum, window=window)
        '''
        
        raise a

        count += 1
        if count > 10:
            raise a        
        
