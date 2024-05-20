import re
import numpy as np
from hashlib import md5
from specutils.io.registers import data_loader
from specutils.io.parsing_utils import read_fileobj_or_hdulist

from Grok.specutils.spectrum import (Spectrum, SpectrumCollection)
from Grok.specutils.utils import concatenate_wat_headers, compute_dispersion, e_flux_to_ivar, get_meta_dict


# Identifiers

def identify_multispec(origin, *args, **kwargs):
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        ctype1 = hdulist[0].header["CTYPE1"].lower()
        wat0_001 = hdulist[0].header["WAT0_001"].lower()
        return (ctype1.startswith("multispe") or wat0_001 == "system=multispec")
    
def identify_pepsi(origin, *args, **kwargs):
    with read_fileobj_or_hdulist(*args, memmap=False, **kwargs) as hdulist:
        return hdulist[0].header.get("INSTRUME") == "PEPSI"


@data_loader(
    "PEPSI", 
    dtype=Spectrum,
    identifier=identify_pepsi, 
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
                "0da149208a3c8ba608226544605ed600": (1, 2, e_flux_to_ivar), # CarPy MIKE product
                "e802331006006930ee0e60c7fbc66cec": (1, 2, e_flux_to_ivar), # Old CarPy MIKE product
                "6b2c2ec1c4e1b122ccab15eb9bd305bc": (1, 2, e_flux_to_ivar), # CarPy MAGE product
                "a4d8f6f51a7260fce1642f7b42012969": (0, 2, e_flux_to_ivar), # IRAF 3 band product
                "148aa0c459c8085f7461a519b1a060e5": (0, None, lambda x: x), # IRAF 1 band product
                "2ab648afed96dcff5ccd10e5b45730c1": (1, 2, e_flux_to_ivar), # DuPont product
            }
            
            try:
                default_bandid_flux, default_bandid_noise, transform = exts[md5_hash]
            except KeyError:
                raise KeyError("Unrecognised multispec type. Cannot identify flux/noise bands. Use `bandid_flux` and `bandid_noise` keywords (1-indexed).")

            bandid_flux = bandid_flux or default_bandid_flux
            bandid_noise = bandid_noise or default_bandid_noise            
        else:
            transform = lambda x: x
            
        flux = hdulist[0].data[bandid_flux - 1]
        if bandid_noise is None:
            ivar = np.ones_like(flux)
        else:
            ivar = transform(hdulist[0].data[bandid_noise - 1])
        
        meta = get_meta_dict(hdulist)
            
    # Ensure λ maps from blue to red direction.
    if np.min(λ[0]) > np.min(λ[-1]):
        λ = λ[::-1]
        flux = flux[:, ::-1]
        ivar = ivar[:, ::-1]
    
    return SpectrumCollection(
        λ=λ,
        flux=flux,
        ivar=ivar,
        meta=meta,
        vacuum=kwargs.get("vacuum", False)
    )    
