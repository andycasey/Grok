
import numpy as np
from typing import Sequence, Optional, Union
from astropy.nddata import NDIOMixin
from astropy.constants import c
from astropy import units as u
from specutils.utils.wcs_utils import air_to_vac

C_KM_S = c.to("km/s").value

# TODO: put elsewhere
def apply_relativistic_doppler_shift(λ, v):
    beta = v / C_KM_S
    return λ * np.sqrt((1 + beta) / (1 - beta))

class Spectrum(NDIOMixin):
    
    def __init__(
        self,
        λ: Sequence[float],
        flux: Sequence[float],
        ivar: Optional[Sequence[float]] = None,
        pixel_flags: Optional[Sequence[int]] = None,
        mask: Optional[Sequence[bool]] = None,
        meta: Optional[dict] = None,
        v_rel: float = 0.0,
        vacuum: Optional[bool] = True,
        continuum: Optional[Union[float, Sequence[float]]] = 1 # should this allow a continuum model of some sort?
    ):
        if ivar is None:
            ivar = np.ones_like(flux)        
        meta = {} if meta is None else meta        
        self.λ = np.array(λ)
        self.flux = np.array(flux)
        self.ivar = np.array(ivar)
        self.mask = np.zeros_like(self.flux, dtype=bool) if mask is None else np.array(mask, dtype=bool)
        self.pixel_flags = pixel_flags
        bad_pixel = ~np.isfinite(self.flux) | ~np.isfinite(self.ivar)
        self.ivar[bad_pixel] = 0        
        self.meta = meta
        self.vacuum = vacuum
        self.v_rel = v_rel
        self.continuum = continuum
        return None
    
            
    @property
    def _relativistic_doppler_shift(self):
        beta = (self.v_rel / C_KM_S)
        return np.sqrt((1 + beta) / (1 - beta))
    
    @property
    def λ_rest(self):
        return self.λ * self._relativistic_doppler_shift
    
    @property
    def λ_rest_vacuum(self):
        return self.λ_vacuum * self._relativistic_doppler_shift
        
    @property
    def λ_vacuum(self):
        if self.vacuum:
            return self.λ
        else:            
            return air_to_vac(self.λ << u.Angstrom).value

    @property
    def rectified_flux(self):
        return self.flux / self.continuum

    @property
    def rectified_ivar(self):
        return self.ivar * self.continuum**2
    
    def __len__(self):
        return self.λ.size

    def apply_velocity_shift(self, v):
        beta = v / C_KM_S
        scale = np.sqrt((1 + beta) / (1 - beta))     
        self.λ *= scale

class SpectrumCollection(Spectrum):
    pass
