import os
import re
import numpy as np
from hashlib import md5
from specutils.io.registers import data_loader
from specutils.io.parsing_utils import read_fileobj_or_hdulist

from Grok.spectrum.spectrum import (Spectrum, SpectrumCollection)
from Grok.spectrum.utils import (
    concatenate_wat_headers, compute_linear_dispersion, compute_dispersion, sigma_to_ivar, get_meta_dict
)

"""
A helpful note for the programmer:

The naming of `Spectrum` and `SpectrumCollection` loosely follows the convention set by astropy.

- Use `Spectrum` for a single one-dimensional spectrum.
- Use `SpectrumCollection` for a collection of one-dimensional spectra that have the same flux shape
  (e.g., echelle orders are often the same shape because they have the same pixels per echelle order,
  even if the wavelength per order is different).

- Use `SpectrumList` for a collection of one-dimensional spectra that have different flux shapes.

** `SpectrumList` is not yet implemented yet **
"""


# Identifiers
def path_given_and_does_not_match(args, pattern):    
    path = args[0]
    return (
        isinstance(path, (str, os.PathLike))
        and 
        (re.match(pattern, os.path.basename(path)) is None)
    )

def _get_instrument(*args, **kwargs):
    with read_fileobj_or_hdulist(*args, **kwargs) as hdulist:
        foo = hdulist[0].header.get("INSTRUME", None)
        print(f"foo = {foo}")
        return foo


def identify_multispec(origin, *args, **kwargs):
    if path_given_and_does_not_match(args, ".+\.fits$"):
        return False
    
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        ctype1 = hdulist[0].header["CTYPE1"].lower()
        wat0_001 = hdulist[0].header["WAT0_001"].lower()
        return (ctype1.startswith("multispe") or wat0_001 == "system=multispec")
    
def identify_pepsi_spectrum(origin, *args, **kwargs):
    return _get_instrument(*args, **kwargs) == "PEPSI"

def identify_generic_spectrum(origin, *args, **kwargs):
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        # Must have CRRVAL1 and CDELT1
        return (
            hdulist[0].header.get("CRVAL1", None) is not None
        and hdulist[0].header.get("CDELT1", None) is not None
        and hdulist[0].data is not None
        )

# This code snippet is a part of the `identify_harps_spectrum` function.
def identify_harps_spectrum(origin, *args, **kwargs):
    if path_given_and_does_not_match(args, "^ADP"):
        return False
        
    return _get_instrument(*args, **kwargs) == "HARPS"

def identify_neid_spectrum(origin, *args, **kwargs):
    if path_given_and_does_not_match(args, ".+\.fits$"):
        return False        
    return _get_instrument(*args, **kwargs) == "NEID"

# Loaders

@data_loader(
    "generic",
    dtype=Spectrum,
    identifier=identify_generic_spectrum,
    extensions=["fits"]
)
def generic(path, **kwargs):
    with read_fileobj_or_hdulist(path, **kwargs) as hdulist:
        f = hdulist[0].data
        ivar = np.ones_like(f)
        mask = np.zeros_like(f, dtype=bool)
        meta = get_meta_dict(hdulist)
        
        λ = compute_linear_dispersion(
            hdulist[0].header["CRVAL1"],
            hdulist[0].header["CDELT1"],
            hdulist[0].header.get("NAXIS1", f.size),
            crpix=hdulist[0].header.get("CRPIX1", 1),
            ltv=hdulist[0].header.get("LTV1", 0)
        )            
            
    return Spectrum(
        λ=λ,
        flux=f,
        ivar=ivar,
        mask=mask,
        meta=meta,
        **kwargs,
    )        


@data_loader(
    "PEPSI", 
    dtype=Spectrum,
    identifier=identify_pepsi_spectrum, 
    extensions=["nor", "avr"]
)
def pepsi(path, **kwargs):
    with read_fileobj_or_hdulist(path, **kwargs) as hdulist:
        λ = hdulist[1].data["Arg"]
        f = hdulist[1].data["Fun"]
        var_f = hdulist[1].data["Var"]
        mask = hdulist[1].data["Mask"]
        meta = get_meta_dict(hdulist)

    # PEPSI spectra are placed in the stellar rest frame by the pipeline
    #v_rel = meta.get("SSBVEL", 0) / 1000.0 # cm/s to km/s
    
    return Spectrum(
        λ=λ,
        flux=f,
        ivar=1/var_f,
        mask=mask,
        meta=meta,
        vacuum=kwargs.get("vacuum", False)
    )

    
@data_loader(
    "multispec",
    dtype=SpectrumCollection,
    identifier=identify_multispec,
    extensions=["fits"]
)
def multispec(path, **kwargs):
    
    # The multispec format fits uses 68, but some files are broken.
    wat_length = kwargs.pop("wat_length", 68)
    
    bandid_flux = kwargs.pop("bandid_flux", None)
    bandid_noise = kwargs.pop("bandid_noise", None)
    
    with read_fileobj_or_hdulist(path, **kwargs) as hdulist:
        
        wat = concatenate_wat_headers(hdulist[0].header, wat_length)

        # Split the concatenated header into individual orders.
        order_mapping = np.array([map(float, each.rstrip('" ').split()) \
                for each in re.split('spec[0-9]+ ?= ?"', wat)[1:]])

        # Parse the order mapping into λ values.
        # Do it this way to ensure ragged arrays work
        num_pixels, num_orders = hdulist[0].header["NAXIS1"], hdulist[0].header["NAXIS2"]
        λ = np.zeros((num_orders, num_pixels), dtype=float) + np.nan
        for j in range(num_orders):
            _ = compute_dispersion(*order_mapping[j])
            λ[j,0:len(_)] = _
            
        # Get the correct extensions.
        if bandid_flux is None or bandid_noise is None:                                        
            md5_hash = md5(";".join([v for k, v in hdulist[0].header.items() \
                                        if k.startswith("BANDID")]).encode("utf-8")).hexdigest()            
            exts = {
                "0da149208a3c8ba608226544605ed600": (1, 2, sigma_to_ivar), # CarPy MIKE product
                "e802331006006930ee0e60c7fbc66cec": (1, 2, sigma_to_ivar), # Old CarPy MIKE product
                "6b2c2ec1c4e1b122ccab15eb9bd305bc": (1, 2, sigma_to_ivar), # CarPy MAGE product
                "a4d8f6f51a7260fce1642f7b42012969": (0, 2, sigma_to_ivar), # IRAF 3 band product
                "148aa0c459c8085f7461a519b1a060e5": (0, None, lambda x: x), # IRAF 1 band product
                "2ab648afed96dcff5ccd10e5b45730c1": (1, 2, sigma_to_ivar), # DuPont product
            }
            
            try:
                default_bandid_flux, default_bandid_noise, transform = exts[md5_hash]
            except KeyError:
                raise KeyError("Unrecognised multispec type. Cannot identify flux/noise bands. Use `bandid_flux` and `bandid_noise` keywords (1-indexed).")

            bandid_flux = bandid_flux or default_bandid_flux
            bandid_noise = bandid_noise or default_bandid_noise            
        else:
            transform = lambda x: x
            
        flux = hdulist[0].data[bandid_flux] 
        if bandid_noise is None:
            ivar = np.ones_like(flux)
        else:
            ivar = transform(hdulist[0].data[bandid_noise])
        
        meta = get_meta_dict(hdulist)
            
    # Ensure λ maps from blue to red direction.
    if np.min(λ[0]) > np.min(λ[-1]):
        ordered = slice(None, None, -1)
        λ = λ[ordered]
        flux = flux[ordered]
        ivar = ivar[ordered]
    
    return SpectrumCollection(
        λ=λ,
        flux=flux,
        ivar=ivar,
        meta=meta,
        vacuum=kwargs.get("vacuum", False)
    )    


@data_loader(
    "HARPS",
    dtype=Spectrum,
    identifier=identify_harps_spectrum,
    extensions=["fits", "fits.gz"]
)
def harps(path, **kwargs):

    with read_fileobj_or_hdulist(path, **kwargs) as hdulist:        
        λ = hdulist[1].data["WAVE"][0]
        flux = hdulist[1].data["FLUX"][0]
        e_flux = hdulist[1].data["ERR"][0]
        if not np.any(np.isfinite(e_flux)):        
            ivar = np.ones_like(flux)
        else:
            ivar = sigma_to_ivar(e_flux)
        
        meta = get_meta_dict(hdulist)
    
    return Spectrum(
        λ=λ,
        flux=flux,
        ivar=ivar,
        meta=meta,
        vacuum=kwargs.get("vacuum", False)
    )        

@data_loader(
    "NEID",
    dtype=SpectrumCollection,
    identifier=identify_neid_spectrum,
    extensions=["fits"]
)
def neid(path, **kwargs):
    with read_fileobj_or_hdulist(path, **kwargs) as hdulist:    
        return SpectrumCollection(
            λ=hdulist[7].data,
            flux=hdulist[1].data,
            ivar=1.0/hdulist[4].data,
            meta=get_meta_dict(hdulist),
            vacuum=True, # THANK FUCK: NEID spectra are in vacuum wavelengths
        )