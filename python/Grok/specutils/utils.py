#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" General utilities for dealing with spectra. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["calculate_fractional_overlap", "find_overlaps"]

import numpy as np
import time

def e_flux_to_ivar(e_flux):
    with np.errstate(divide='ignore'):
        return e_flux**-2

def get_meta_dict(hdulist):
    return dict(hdulist[0].header)


def compute_dispersion(aperture, beam, dispersion_type, dispersion_start,
    mean_dispersion_delta, num_pixels, redshift, aperture_low, aperture_high,
    weight=1, offset=0, function_type=None, order=None, Pmin=None, Pmax=None,
    *coefficients):
    """
    Compute a dispersion mapping from a IRAF multi-spec description.

    :param aperture:
        The aperture number.

    :param beam:
        The beam number.

    :param dispersion_type:
        An integer representing the dispersion type:

        0: linear dispersion
        1: log-linear dispersion
        2: non-linear dispersion

    :param dispersion_start:
        The value of the dispersion at the first physical pixel.

    :param mean_dispersion_delta:
        The mean difference between dispersion pixels.

    :param num_pixels:
        The number of pixels.

    :param redshift:
        The redshift of the object. This is accounted for by adjusting the
        dispersion scale without rebinning:

        >> dispersion_adjusted = dispersion / (1 + redshift)

    :param aperture_low:
        The lower limit of the spatial axis used to compute the dispersion.

    :param aperture_high:
        The upper limit of the spatial axis used to compute the dispersion.

    :param weight: [optional]
        A multiplier to apply to all dispersion values.

    :param offset: [optional]
        A zero-point offset to be applied to all the dispersion values.

    :param function_type: [optional]
        An integer representing the function type to use when a non-linear 
        dispersion mapping (i.e. `dispersion_type = 2`) has been specified:

        1: Chebyshev polynomial
        2: Legendre polynomial
        3: Cubic spline
        4: Linear spline
        5: Pixel coordinate array
        6: Sampled coordinate array

    :param order: [optional]
        The order of the Legendre or Chebyshev function supplied.

    :param Pmin: [optional]
        The minimum pixel value, or lower limit of the range of physical pixel
        coordinates.

    :param Pmax: [optional]
        The maximum pixel value, or upper limit of the range of physical pixel
        coordinates.

    :param coefficients: [optional]
        The `order` number of coefficients that define the Legendre or Chebyshev
        polynomial functions.

    :returns:
        An array containing the computed dispersion values.
    """

    if dispersion_type in (0, 1):
        # Simple linear or logarithmic spacing
        dispersion = \
            dispersion_start + np.arange(num_pixels) * mean_dispersion_delta

        if dispersion_start == 1:
            dispersion = 10.**dispersion

    elif dispersion_type == 2:
        # Non-linear mapping.
        if function_type is None:
            raise ValueError("function type required for non-linear mapping")
        elif function_type not in range(1, 7):
            raise ValueError(
                "function type {0} not recognised".format(function_type))

        assert Pmax == int(Pmax), Pmax; Pmax = int(Pmax)
        assert Pmin == int(Pmin), Pmin; Pmin = int(Pmin)

        if function_type == 1:
            # Chebyshev polynomial.
            if None in (order, Pmin, Pmax, coefficients):
                raise TypeError("order, Pmin, Pmax and coefficients required "
                                "for a Chebyshev or Legendre polynomial")

            order = int(order)
            n = np.linspace(-1, 1, Pmax - Pmin + 1)
            temp = np.zeros((Pmax - Pmin + 1, order), dtype=float)
            temp[:, 0] = 1
            temp[:, 1] = n
            for i in range(2, order):
                temp[:, i] = 2 * n * temp[:, i-1] - temp[:, i-2]
            
            for i in range(0, order):
                temp[:, i] *= coefficients[i]

            dispersion = temp.sum(axis=1)


        elif function_type == 2:
            # Legendre polynomial.
            if None in (order, Pmin, Pmax, coefficients):
                raise TypeError("order, Pmin, Pmax and coefficients required "
                                "for a Chebyshev or Legendre polynomial")

            Pmean = (Pmax + Pmin)/2
            Pptp = Pmax - Pmin
            x = (np.arange(num_pixels) + 1 - Pmean)/(Pptp/2)
            p0 = np.ones(num_pixels)
            p1 = mean_dispersion_delta

            dispersion = coefficients[0] * p0 + coefficients[1] * p1
            for i in range(2, int(order)):
                if function_type == 1:
                    # Chebyshev
                    p2 = 2 * x * p1 - p0
                else:
                    # Legendre
                    p2 = ((2*i - 1)*x*p1 - (i - 1)*p0) / i

                dispersion += p2 * coefficients[i]
                p0, p1 = (p1, p2)

        elif function_type == 3:
            # Cubic spline.
            if None in (order, Pmin, Pmax, coefficients):
                raise TypeError("order, Pmin, Pmax and coefficients required "
                                "for a cubic spline mapping")
            s = (np.arange(num_pixels, dtype=float) + 1 - Pmin)/(Pmax - Pmin) \
              * order
            j = s.astype(int).clip(0, order - 1)
            a, b = (j + 1 - s, s - j)
            x = np.array([
                a**3,
                1 + 3*a*(1 + a*b),
                1 + 3*b*(1 + a*b),
                b**3])
            dispersion = np.dot(np.array(coefficients), x.T)

        else:
            raise NotImplementedError("function type not implemented yet")

    else:
        raise ValueError(
            "dispersion type {0} not recognised".format(dispersion_type))

    # Apply redshift correction.
    dispersion = weight * (dispersion + offset) / (1 + redshift)
    return dispersion

def concatenate_wat_headers(header, wat_length=68, wat_prefix="WAT2_"):
    i, wat, key_fmt = (1, str(""), "{wat_prefix}{i:03d}")
    while True:
        value = header.get(key_fmt.format(wat_prefix=wat_prefix, i=i))
        if value is None: 
            return wat
        wat += value + (" "  * (wat_length - len(value)))
        i += 1


def calculate_fractional_overlap(interest_range, comparison_range):
    """
    Calculate how much of the range of interest overlaps with the comparison
    range.
    """

    if not (interest_range[-1] >= comparison_range[0] \
        and comparison_range[-1] >= interest_range[0]):
        return 0.0 # No overlap

    elif   (interest_range[0] >= comparison_range[0] \
        and interest_range[-1] <= comparison_range[-1]):
        return 1.0 # Total overlap 

    else:
        # Some overlap. Which side?
        if interest_range[0] < comparison_range[0]:
            # Left hand side
            width = interest_range[-1] - comparison_range[0]

        else:
            # Right hand side
            width = comparison_range[-1] - interest_range[0]
        return width/np.ptp(interest_range) # Fractional overlap


def find_overlaps(spectra, dispersion_range, return_indices=False):
    """
    Find spectra that overlap with the dispersion range given. Spectra are
    returned in order of how much they overlap with the dispersion range.

    :param spectra:
        A list of spectra.

    :param dispersion_range:
        A two-length tuple containing the start and end wavelength.

    :param return_indices: [optional]
        In addition to the overlapping spectra, return their corresponding
        indices.

    :returns:
        The spectra that overlap with the dispersion range provided, ranked by
        how much they overlap. Optionally, also return the indices of those
        spectra.
    """

    fractions = np.array([
        calculate_fractional_overlap(s.dispersion, dispersion_range) \
            for s in spectra])

    N = (fractions > 0).sum()
    indices = np.argsort(fractions)[::-1]
    overlaps = [spectra[index] for index in indices[:N]]

    """
    A faster, alternative method if sorting is not important:
    # http://stackoverflow.com/questions/325933/determine-whether-two-date-ranges-overlap/325964#325964    
    overlaps, indices = zip(*[(spectrum, i) \
        for i, spectrum in enumerate(spectra) \
            if  spectrum.dispersion[-1] >= dispersion_range[0] \
            and dispersion_range[-1] >= spectrum.dispersion[0]])
    """

    return overlaps if not return_indices else (overlaps, indices[:N])

def calculate_noise(y, window, method='madstd', verbose=True, full_output=False):
    """
    Calculate SNR a posteriori.
    Very inefficient but who cares!
    :param y:
        array of flux values to calculate SNR.
        Assumes equal or slowly changing x sampling
    :param window:
        pixel window to use
    :param method:
        madstd (default): use astropy.stats.mean_absolute_deviation
        betasigma: use PyAstronomy BSEqSamp
    :param verbose:
        if True (default), prints how long it takes to calculate the SNR.
    """
    assert method in ['madstd', 'betasigma']

    from scipy import signal
    ## Signal
    S = signal.medfilt(y, window)
    
    ## Noise Estimator
    if method == 'betasigma':
        from PyAstronomy import pyasl
        beq = pyasl.BSEqSamp()
        N = 1; j = 1
        def noise_estimator(x):
            return beq.betaSigma(x, N, j, returnMAD=True)[0]
    elif method == 'madstd':
        try:
            from astropy.stats import mean_absolute_deviation
        except ImportError:
            def noise_estimator(x):
                return np.nanmedian(np.abs(x-np.nanmedian(x)))
        else:
            noise_estimator = lambda x: mean_absolute_deviation(x, ignore_nan=True)
    
    ## Calculate noise
    noise = np.zeros(len(y))
    if verbose: start = time.time()
    # middle: calculate rolling value
    for i in range(window, len(noise)-window):
        noise[i] = noise_estimator(y[i-window:i+window])
    # edges: just use the same noise value
    noise[:window] = noise_estimator(y[:window])
    noise[-window:] = noise_estimator(y[-window:])
    
    if verbose:
        print("Noise estimate of {} points with window {} took {:.1f}s".format(len(noise), window, time.time()-start))

    # https://en.wikipedia.org/wiki/Median_absolute_deviation
    # sigma = 1.4826 * MAD
    noise *= 1.4826
    
    if full_output:
        SNR = S/noise
        return noise, S, SNR
    return noise

def calculate_snr(*args, **kwargs):
    kwargs = kwargs.copy()
    kwargs["full_output"] = True
    noise, S, SNR = calculate_noise(*args, **kwargs)
    return SNR
