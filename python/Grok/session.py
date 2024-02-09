import os
import numpy as np
from typing import Union, Sequence, Type, Optional
from functools import cached_property


from Grok.specutils import Spectrum1D
from Grok.synthesis import BaseKorg, Korg

# Move to utilities:
        
def expand_path(path):
    return os.path.abspath(os.path.expanduser(path))

def expand_input_paths(input_paths):
    if isinstance(input_paths, (str, bytes)):
        input_paths = [input_paths]
    return tuple(map(expand_path, input_paths))


def get_closest_spectrum_index(session, wavelength):
    """Return the spectrum in the session closest to the wavelength."""
    mean_wavelengths = [np.mean(wl) for wl, *_ in session.spectra]
    index = np.argmin(np.abs(np.array(mean_wavelengths) - wavelength))
    return index




class Session:
    
    def __init__(
        self,
        input_paths: Sequence[Union[str, bytes, os.PathLike]],
        synthesis: Optional[Type[BaseKorg]] = Korg
    ):
        self.input_paths = expand_input_paths(input_paths)
        self._synthesis = synthesis
        return None
    
    @cached_property
    def synthesis(self):
        # If ._synthesis is a class, instantiate it
        if isinstance(self._synthesis, type): 
            return self._synthesis()
        else:
            return self._synthesis
    
    
    @cached_property
    def spectra(self):
        spectra = []
        for input_path in self.input_paths:            
            try:
                wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec(input_path)
            except:
                #    wavelength, flux, ivar, meta = Spectrum1D.read_apogee(input_path)
                from astropy.io import fits
                with fits.open(input_path) as image:
                    flux = image[0].data
                    wavelength = image[0].header["CRVAL1"] + image[0].header["CDELT1"] * np.arange(image[0].header["NAXIS1"])
                    ivar = np.ones_like(flux)
                    
                wavelength = [wavelength]            
                flux = [flux]
                ivar = [ivar]
                meta = {"input_path": input_path}
                
            for i in range(len(wavelength)):
                spectra.append([wavelength[i], flux[i], ivar[i], meta])
        return spectra
    
    
    def _get_closest_spectrum_index(self, wavelength):
        return get_closest_spectrum_index(self, wavelength)
        
    
    
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
            
    
    #def measure_relative_velocity(self, index, template_path, continuum_kwargs=None, v_lim=(-250, 250), deg=1, N=10, **kwargs):
        
        

    # Curve-of-growth analysis
    def cog_read_linelist(self, line_list_path, format="vald", **kwargs):            
        self.cog_linelist = self.synthesis.read_linelist(
            line_list_path, 
            format=format,
            **kwargs
        )
        return self.cog_linelist
    
    
    def cog_fit_profiles(self):
        # just assign EWs to everything in the .linelist        
        ...
        
    
    def cog_interpolate_atmosphere(self, Teff, logg, metals_H, alpha_H):
        self.cog_A_X = self.synthesis.format_A_X(metals_H, alpha_H)

        self.cog_atm = self.synthesis.interpolate_marcs(
            Teff, logg, self.cog_A_X
        )
        
        
    def cog_ews_to_abundances(self, ews, **kwargs):
        return self.synthesis.ews_to_abundances(
            self.cog_atm,
            self.cog_linelist,
            self.cog_A_X,
            ews,
            **kwargs
        )
        
    
    
    def cog_solve_stellar_parameters(
        self,
        ews,
        Teff0=5000.0, 
        logg0=3.5, 
        vmic0=1.0, 
        metallicity0=0.0, 
        **kwargs
    ):
        return self.synthesis.ews_to_stellar_parameters(
            self.cog_linelist, 
            ews, 
            Teff0=Teff0,
            logg0=logg0,
            vmic0=vmic0,
            metallicity0=metallicity0,
            **kwargs
        )
