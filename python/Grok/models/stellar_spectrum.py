
import numpy as np
from scipy import optimize as op
from itertools import accumulate, cycle
from types import MappingProxyType
from typing import Sequence, Dict, Union, Optional, Tuple

from Grok.models.basis import LinearBasis
from Grok.spectrum.utils import ivar_to_sigma
from Grok.utils import overlap

from Grok.spectrum import Spectrum, SpectrumCollection

class Model:
    ...
    

class LinearStellarSpectrumModel(Model):
    
    def __init__(
        self,
        stellar_basis: LinearBasis, 
        telluric_basis: Optional[LinearBasis] = None,
        continuum_basis: Optional[LinearBasis] = None,
        mask_regions: Optional[Sequence[Tuple[float, float]]] = None,
    ):
        """
        Model a stellar spectrum as a (nearly) linear combination of basis spectra.
        
        :param stellar_basis:
            The basis for the stellar spectrum.
        
        :param telluric_basis: [optional]
            The basis for the telluric spectrum.
            
        :param continuum_basis: [optional]
            The basis for the continuum.
        """
        
        self.bases = MappingProxyType(dict([
            (k, v) for k, v in [
                ("stellar", stellar_basis),
                ("telluric", telluric_basis),
                ("continuum", continuum_basis)
            ] if v is not None        
        ]))
        self.mask_regions = mask_regions
        return None
    
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
            
    
    def prepare_spectra(self, spectra):
        if isinstance(spectra, Spectrum):
            spectra = [spectra]

        S = len([s.flux.shape[0] if isinstance(s, SpectrumCollection) else 1 for s in spectra])
        P = sum(tuple(map(len, spectra)))
        
        λ, flux, ivar, oi = (np.empty(P), np.empty(P), np.empty(P), np.empty(P, dtype=int))

        i, si = (0, 0)
        for spectrum in spectra:
            if isinstance(spectrum, SpectrumCollection):
                for j in range(spectrum.flux.shape[0]):
                    sliced = slice(si, si + spectrum.λ[j].size)
                    oi[sliced] = i
                    λ[sliced] = spectrum.λ_rest_vacuum[j]
                    flux[sliced] = spectrum.flux[j]
                    ivar[sliced] = spectrum.ivar[j]
                    si += spectrum.λ[j].size
                    i += 1
            else:
                sliced = slice(si, si + spectrum.λ.size)
                oi[sliced] = i
                λ[sliced] = spectrum.λ_rest_vacuum
                flux[sliced] = spectrum.flux
                ivar[sliced] = spectrum.ivar
                si += spectrum.λ.size
                i += 1                

        pi = np.argsort(λ) # pixel indices
        λ, flux, ivar, oi = (λ[pi], flux[pi], ivar[pi], oi[pi])
        
        zero_out_bad_pixels(flux, ivar)
        mask = ~self.get_mask(λ)
        
        return (λ, flux, ivar, mask, oi, pi, S)
    
    
    def prepare_data_arrays(self, λ, flux, ivar: Optional[Sequence[float]]):
        λ = np.atleast_1d(λ).flatten()
        flux = np.atleast_1d(flux).flatten()
        if ivar is None:
            ivar = np.ones_like(flux)
        else:
            ivar = np.atleast_1d(ivar).flatten()
        oi = np.zeros_like(λ, dtype=int)
        pi = np.argsort(λ)
        λ, flux, ivar, oi = (λ[pi], flux[pi], ivar[pi], oi[pi])
        zero_out_bad_pixels(flux, ivar)        
        P = λ.size        
        mask = ~self.get_mask(λ)
        return (λ, flux, ivar, mask, oi, pi, 1)        
    
    
    def _fit(
        self, 
        λ: Sequence[float], 
        flux: Sequence[float],
        ivar: Optional[Sequence[float]] = None,
        x0: Optional[Sequence[float]] = None,
        **kwargs
    ):
        """
        Fit the model to the given data, assuming it is a single spectrum order.
        
        :param λ:
            The rest vacuum wavelengths of the spectrum.
        
        :param flux:
            The flux of the spectrum.
        
        :param ivar: [optional]
            The inverse variance of the spectrum.
            
        :param x0: [optional]
            The initial guess for the model parameters.
        """

        (λ, flux, ivar, mask, oi, pi, S) = (*args, pi, S) = self.prepare_data_arrays(λ, flux, ivar)
        
        f, g, p0 = self.get_forward_model_and_initial_guess(*args, x0)

        self.θ, self.Σ = op.curve_fit(
            f,
            λ[mask],
            flux[mask],
            p0=p0,
            jac=g,
            sigma=ivar_to_sigma(ivar[mask]),
            bounds=self.get_bounds(oi)
        )
        y, rectified_stellar_flux, rectified_telluric_flux, continuum = self(λ, oi, *self.θ, full_output=True)                
        
        return (λ, oi, y, rectified_stellar_flux, rectified_telluric_flux, continuum)

    

    def get_bounds(self, oi):
        """
        Get the bounds for the model parameters.
        
        :param oi:
            The order indices: the index of the order in the spectrum.
        """
        
        bounds = list(self.bases["stellar"].get_bounds())
        
        if "telluric" in self.bases:
            bounds.extend(self.bases["telluric"].get_bounds())
        
        if "continuum" in self.bases:         
            continuum_basis = self.bases["continuum"]
            if not isinstance(continuum_basis, (list, tuple)):
                continuum_basis = [continuum_basis]
            O = 1 + np.max(oi)
            for o, base in zip(range(O), cycle(continuum_basis)):
                bounds.extend(base.get_bounds())
        
        return np.array(bounds).T
    
    
    def get_forward_model_and_initial_guess(self, λ, flux, ivar, mask, oi, p0: Optional[Sequence[float]] = None):
        """
        Get the forward model and initial guess for the model parameters.
        
        :param λ:
            The rest vacuum wavelengths of the spectrum.
        
        :param flux:
            The flux of the spectrum.
        
        :param ivar:
            The inverse variance of the spectrum.
        
        :param mask:
            The mask for the spectrum.
        
        :param oi:
            A sequence of order indices: the index of the order in the spectrum.
        
        :param p0: [optional]
            The initial guess for the model parameters.
            
        :returns:
            A tuple of (f, g, p0), where f is the forward model, g is the Jacobian of the forward model, and p0
            is the initial guess for the model parameters.
        """
        
        x0 = list(self.bases["stellar"].get_initial_guess(λ[mask], flux[mask], ivar[mask]))    
        A_stellar = self.bases["stellar"].get_design_matrix(λ[mask])
        
        # Tellurics
        try:
            telluric_basis = self.bases["telluric"]
        except KeyError:
            _A = -A_stellar
        else:
            A_telluric = telluric_basis.get_design_matrix(λ[mask])
            _A = -np.hstack([A_stellar, A_telluric])
            x0.extend(telluric_basis.get_initial_guess(λ[mask], flux[mask], ivar[mask]))

        f, g = (None, None)
            
        # Continuum
        try:
            continuum_basis = self.bases["continuum"]
        except KeyError:
            # No continuum modelling.
            
            def f(λ, *θ, full_output=False):
                if not full_output:
                    return 1 + _A @ θ
                else:
                    rectified_stellar_flux = 1 + _A[:, :A_stellar.shape[1]] @ θ[:A_stellar.shape[1]]
                    if "telluric" in self.bases:
                        rectified_telluric_flux = 1 + _A[:, A_stellar.shape[1]:] @ θ[A_stellar.shape[1]:]
                    else:
                        rectified_telluric_flux = 1                    
                    return (1 + _A @ θ, rectified_stellar_flux, rectified_telluric_flux, 1)
            
            def g(λ, *θ):
                return _A

            if p0 is None:
                p0 = x0 
            return (f, g, p0)                            
        
        else:
            O = 1 + int(np.max(oi))
            if isinstance(continuum_basis, (list, tuple)):
                if len(continuum_basis) != O:
                    raise ValueError(f"Continuum basis must be a list of length O ({len(continuum_basis)} != {O})")
                n_continuum_parameters = sum(b.n_parameters for b in continuum_basis)
            else:
                n_continuum_parameters = continuum_basis.n_parameters * O
                continuum_basis = [continuum_basis] * O
                                
            A_continuum = np.zeros((λ.size, n_continuum_parameters))
            si, order_masks = (0, [])
            for o, base in enumerate(continuum_basis):
                no = base.n_parameters
                order_mask = (oi == o) * mask
                order_masks.append(order_mask)
                A_order = base.get_design_matrix(λ[order_mask])
                A_continuum[order_mask, si:si+no] = A_order
                si += no
                if p0 is None:
                    x0.extend(base.get_initial_guess(λ[order_mask], flux[order_mask], ivar[order_mask], A=A_order))
            
            n = _A.shape[1]
                                    
            def f(λ, *θ, full_output=False):
                y = (1 + _A @ θ[:n]) * (A_continuum @ θ[n:])                      
                if not full_output:
                    return y
                else:
                    rectified_stellar_flux = 1 + _A[:, :A_stellar.shape[1]] @ θ[:A_stellar.shape[1]]
                    if "telluric" in self.bases:
                        rectified_telluric_flux = 1 + _A[:, A_stellar.shape[1]:] @ θ[A_stellar.shape[1]:A_stellar.shape[1] + A_telluric.shape[1]]
                    else:
                        rectified_telluric_flux = 1
                    continuum = A_continuum @ θ[-A_continuum.shape[1]:]
                    
                    return (y, rectified_stellar_flux, rectified_telluric_flux, continuum)
            
            def g(λ, *θ):
                return np.hstack([
                    _A * (A_continuum @ θ[n:])[:, np.newaxis],
                    (1 + _A @ θ[:n])[:, np.newaxis] * A_continuum
                ])
                
            if p0 is None:
                p0 = x0 
            return (f, g, p0)

        
    def __call__(self, λ, oi, *θ, **kwargs):        
        dummy = np.ones_like(λ)        
        f, g, _ = self.get_forward_model_and_initial_guess(λ, dummy, dummy, dummy.astype(bool), oi, 1)
        return f(λ, *θ, **kwargs)
        
    
    def fit(self, spectra: Union[Spectrum, Sequence[Spectrum]], x0: Optional[Sequence[float]] = None, **kwargs):
        """
        Fit a sequence of spectra of the same star.
        
        :param spectra:
            The spectra to fit.
        
        :param x0: [optional]
            The initial guess for the model parameters.
        """
        
        (λ, flux, ivar, mask, oi, *_) = (*args, pi, S) = self.prepare_spectra(spectra)                
        f, g, p0 = self.get_forward_model_and_initial_guess(*args, x0)

        self.θ, self.Σ = op.curve_fit(
            f,
            λ[mask],
            flux[mask],
            p0=p0,
            jac=g,
            sigma=ivar_to_sigma(ivar[mask]),
            bounds=self.get_bounds(oi)
        )
                
        y, rectified_stellar_flux, rectified_telluric_flux, continuum = self(λ, oi, *self.θ, full_output=True)                
        
        return (λ, oi, y, rectified_stellar_flux, rectified_telluric_flux, continuum)
        """
        fig, ax = plt.subplots()
        for i in range(1 + np.max(oi)):
            mask = (oi == i)            
            ax.plot(λ[mask], flux[mask], c="k") 
            ax.plot(λ[mask], (rectified_stellar_flux * continuum)[mask], c="r")
            ax.plot(λ[mask], (rectified_telluric_flux * continuum)[mask], c="tab:blue")
        """
        

def zero_out_bad_pixels(flux, ivar):
    bad_pixel = (~np.isfinite(flux)) | (~np.isfinite(ivar)) | (flux <= 0)
    ivar[bad_pixel] = 0
    return None
    
    
        
    
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

        
    from Grok.spectrum import Spectrum        
    from Grok.models.basis import NMFSpectralBasis, PolynomialBasis, FourierBasis
    from Grok.models.profile import VoigtProfile
    from Grok.models.spectral_line import SpectralLineModel
    from Grok.utils import overlap, to_contiguous_regions
    from copy import copy

    
    import h5py as h5
    with h5.File("Grok/bosz-highres-optical-basis-vectors.h5") as fp:
        foo = (fp["vacuum_wavelength"][:], fp["stellar_basis_vectors"][:].T)
        stellar_base = NMFSpectralBasis(fp["vacuum_wavelength"][:], fp["stellar_basis_vectors"][:].T, Ro=27_000, Ri=80_000)
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
        #telluric=telluric_base,
        continuum=FourierBasis(5)
    )
    
    model = LinearStellarSpectrumModel(bases)# mask_regions=mask_regions)
    
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
        
