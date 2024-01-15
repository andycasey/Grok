import numpy as np

def fit_polynomial(wavelength, flux, ivar, mask=None, deg=2):
    if mask is not None:
        raise NotImplementedError("Masking not implemented yet.")
    
    # construct design matrix for fitting polynomial
    A = np.vander(wavelength, deg+1)
    C = np.diag(ivar)
    ATCinvA = np.dot(A.T, np.dot(C, A))
    ATCinvY = np.dot(A.T, np.dot(C, flux))
    X = np.linalg.solve(ATCinvA, ATCinvY)
    
    continuum = np.polyval(X, wavelength)
    return (continuum, X)
    