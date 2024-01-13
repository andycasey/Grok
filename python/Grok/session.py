import os
import numpy as np
from typing import Union, Sequence, Type, Optional
from functools import cached_property


from specutils import Spectrum1D
from synthesis import BaseKorg, Korg

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
        print(self, hash(self), "loading spectra")
        spectra = []
        for input_path in self.input_paths:            
            wavelength, flux, ivar, meta = Spectrum1D.read_fits_multispec(input_path)
            for i in range(len(wavelength)):
                spectra.append([wavelength[i], flux[i], ivar[i], meta])
        return spectra
    
    
    def _get_closest_spectrum_index(self, wavelength):
        return get_closest_spectrum_index(self, wavelength)
        
        

    # Curve-of-growth analysis
    def cog_read_linelist(self, line_list_path, format="vald"):
        self.cog_linelist = self.synthesis.read_linelist(
            line_list_path, 
            format=format,
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
        
        
    def cog_ews_to_abundances(self, ews):
        self.synthesis.ews_to_abundances(
            self.cog_atm,
            self.cog_linelist,
            self.cog_A_X,
            ews
        )
        
    
    
    def cog_solve_stellar_parameters(
        self,
        x0=None,
        callback=None
    ):
        ...
    
    