
import numpy as np
from scipy.signal.windows import tukey
from scipy.optimize import leastsq


def template_logwl_resample(
    wavelength, 
    flux,
    template_wavelength,
    template_flux,
    wblue=None, 
    wred=None,
    delta_log_wavelength=None,
):
    """
    Resample a spectrum and template onto a common log-spaced spectral grid.
    """

    if wblue:
        w0 = np.log10(wblue)
    else:
        ws0 = np.log10(wavelength[0])
        wt0 = np.log10(template_wavelength[0])
        w0 = min(ws0, wt0)

    if wred:
        w1 = np.log10(wred)
    else:
        ws1 = np.log10(wavelength[-1])
        wt1 = np.log10(template_wavelength[-1])
        w1 = max(ws1, wt1)

    if delta_log_wavelength is None:
        ds = np.log10(wavelength[1:]) - np.log10(wavelength[:-1])
        dw = ds[np.argmin(ds)]
    else:
        dw = delta_log_wavelength

    nsamples = int((w1 - w0) / dw)

    log_wave_array = np.ones(nsamples) * w0
    for i in range(nsamples):
        log_wave_array[i] += dw * i

    # Build the corresponding wavelength array
    wave_array = np.power(10., log_wave_array)

    # Resample spectrum and template into wavelength array so built
    resampled_flux = np.interp(wave_array, wavelength, flux, left=0, right=0)
    resampled_template_flux = np.interp(wave_array, template_wavelength, template_flux, left=0, right=0)
    
    return (wave_array, resampled_flux, resampled_template_flux)



def cross_correlate(wavelength, flux, template_wavelength, template_flux, apodization_window=0.2, resample=True):
    if resample:
        if resample is True:
            resample_kwargs = dict()  # use defaults
        else:
            resample_kwargs = resample
        resampled_wavelength, resampled_flux, resampled_template_flux = template_logwl_resample(
            wavelength, flux, template_wavelength, template_flux, **resample_kwargs)

    else:
        resampled_wavelength = wavelength
        resampled_flux = flux
        resampled_template_flux = template_flux

    resampled_flux, resampled_template_flux = _apodize(resampled_flux, resampled_template_flux, apodization_window)

    resampled_flux -= np.nanmean(resampled_flux)
    resampled_template_flux -= np.nanmean(resampled_template_flux)

    corr = np.correlate(resampled_flux, resampled_template_flux, mode='full')

    delta_log_wave = np.log10(resampled_wavelength[1]) - np.log10(resampled_wavelength[0])
    deltas = (np.array(range(len(corr))) - len(corr)/2 + 0.5) * delta_log_wave
    lags = np.power(10., deltas) - 1.
    return (corr, lags * 299792.458)



def measure_relative_velocity(wavelength, flux, template_wavelength, template_flux, v_lim=(-250, 250), deg=1, N=10, **kwargs):
    
    corr, lags = cross_correlate(wavelength, flux, template_wavelength, template_flux, **kwargs)
    
    si, ei = lags.searchsorted(v_lim)
    
    theta = np.polyfit(lags[si:ei], corr[si:ei], deg=deg)
    
    norm_corr = corr - np.polyval(theta, lags)
    
    pi = si + np.argmax(norm_corr[si:ei])
    v_rel = lags[pi]

    # Get initial guess of peak.
    p0 = np.array([v_rel, np.max(norm_corr[si:ei]), 10])# + [0] * deg)

    gaussian = lambda p, x: p[1] * np.exp(-(x - p[0])**2 / (2.0 * p[2]**2))
    errfunc = lambda p, x, y: y - gaussian(p, x)

    # only consider +/- N pixels
    ps, pe = (pi - N, pi + N + 1)
    try:
        p1, ier = leastsq(errfunc, p0.copy(), args=(lags[ps:pe], norm_corr[ps:pe]))
    except:
        raise 
    
    v_rel, e_v_rel = (p1[0], p1[2]) 
    chi2 = -2 * np.max(norm_corr[si:ei])
    
    return (v_rel, e_v_rel, chi2, lags, corr, norm_corr)
    


def _apodize(spectrum, template, apodization_window):
    # Apodization. Must be performed after resampling.
    if apodization_window is None:
        clean_spectrum = spectrum
        clean_template = template
    else:
        if callable(apodization_window):
            window = apodization_window
        else:
            def window(wlen):
                return tukey(wlen, alpha=apodization_window)
        clean_spectrum = spectrum * window(len(spectrum))
        clean_template = template * window(len(template))

    return clean_spectrum, clean_template
