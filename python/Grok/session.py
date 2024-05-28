import os
import numpy as np
from typing import Union, Sequence, Type, Optional
from functools import cached_property

from Grok.utils import expand_path
from Grok.spectrum import Spectrum, SpectrumCollection
from Grok.synthesis import BaseKorg, Korg


class Session:
    
    def __init__(
        self,
        input_paths: Sequence[Union[str, bytes, os.PathLike]],
        synthesis: Optional[Type[BaseKorg]] = Korg
    ):
        # Read spectra.
        self.spectra, self.n_orders_per_spectrum = _read_spectra(input_paths)
        self._synthesis = synthesis
        return None

    
    @cached_property
    def synthesis(self):
        return self._synthesis() if isinstance(self._synthesis, type) else self._synthesis
    
    
    def get_spectrum_indices(self, spectral_order):
        N = np.cumsum(self.n_orders_per_spectrum)
        spectrum_index = N.searchsorted(spectral_order + 1)
        order_index = int(spectral_order - np.sum(self.n_orders_per_spectrum[:spectrum_index]))
        if order_index == 0 and self.n_orders_per_spectrum[spectrum_index] == 1:
            order_index = None
        return (spectrum_index, order_index)
        

    def get_spectral_order(self, spectral_order=0):
        spectrum_index, order_index = self.get_spectrum_indices(spectral_order)
        spectrum = self.spectra[spectrum_index]
        λ, flux, ivar = (spectrum.λ, spectrum.flux, spectrum.ivar)
        if order_index is not None:
            λ, flux, ivar = (λ[order_index], flux[order_index], ivar[order_index])
        return (λ, flux, ivar, spectrum.meta)            
        
    
    def _get_closest_spectrum_index(self, wavelength):
        """Return the spectrum in the session closest to the wavelength."""
        mean_wavelengths = [np.mean(wl) for wl, *_ in self.spectra]
        return np.argmin(np.abs(np.array(mean_wavelengths) - wavelength))

            
    def get_spectrum(self, index=0, rest_frame=False, rectified=False, **kwargs):
        
        _wavelength, _flux, _ivar, meta = self.spectra[index]

        continuum = 1 if not rectified else self.get_continuum(index, **kwargs)            
        v_rel = self.get_v_rel(index) if rest_frame else 0
            
        wavelength = _wavelength * (1 - v_rel / 299792.458)
        flux = _flux / continuum
        ivar = _ivar * continuum**2
        
        return (wavelength, flux, ivar, continuum, v_rel, meta)
    
    
    def get_continuum(self, index=0, **kwargs):
        return 1
    
    

    

    def set_v_rel(self, v_rel):
        self.v_rel = v_rel
        return None
    
    
    def get_v_rel(self, index=None):
        try:
            return self.v_rel[index]
        except:
            try:
                return self.v_rel
            except:
                return 0
    
    # Fit NMF model + continuum?
    
    
    # Fit EW of some line, perhaps given some NMF model for setting masks?
    


def _read_spectra(input_paths):
    spectra, n_orders_per_spectrum = ([], [])
    for path in input_paths:
        for kind in (Spectrum, SpectrumCollection):
            try:
                spectrum = kind.read(expand_path(path))
                n = spectrum.flux.shape[0] if isinstance(spectrum, SpectrumCollection) else 1                
                spectra.append(spectrum)
                n_orders_per_spectrum.append(n)                
            except:
                continue
            else:
                break
        else:
            raise ValueError(f"Could not read {path}")
    return (spectra, tuple(n_orders_per_spectrum))
