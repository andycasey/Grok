
import numpy as np
from typing import Sequence, Optional, Union
from astropy.nddata import NDIOMixin
from astropy.constants import c
from astropy import units as u
from specutils.utils.wcs_utils import air_to_vac

C_KM_S = c.to("km/s").value

MEDIUMS = ("air", "vacuum")

class Spectrum(NDIOMixin):
    
    def __init__(
        self,
        λ: Sequence[float],
        flux: Sequence[float],
        ivar: Optional[Sequence[float]] = None,
        pixel_flags: Optional[Sequence[int]] = None,
        meta: Optional[dict] = None,
        v_rel: float = 0.0,
        medium: str = "vacuum",
        continuum: Optional[Union[float, Sequence[float]]] = 1
    ):
        if ivar is None:
            ivar = np.ones_like(flux)        
        meta = {} if meta is None else meta        
        medium = medium.strip().lower()
        if medium not in ("air", "vacuum"):
            raise ValueError(f"medium must be either 'air' or 'vacuum'")
        self.λ = np.array(λ)
        self.flux = np.array(flux)
        self.ivar = np.array(ivar)
        self.pixel_flags = pixel_flags
        bad_pixel = ~np.isfinite(self.flux) | ~np.isfinite(self.ivar)
        self.ivar[bad_pixel] = 0        
        self.meta = meta
        self.medium = medium
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
        f = air_to_vac if self.medium == "air" else (lambda x: x)        
        return f(self.λ << u.Angstrom).value

    @property
    def rectified_flux(self):
        return self.flux / self.continuum

    @property
    def rectified_ivar(self):
        return self.ivar * self.continuum**2
    
    def __len__(self):
        return self.λ.size



from specutils.io.registers import data_loader


@data_loader("apogee", dtype=Spectrum, verbose=True)
def read_apogee(*args):
    print(f"here we are {args}")
    raise a
