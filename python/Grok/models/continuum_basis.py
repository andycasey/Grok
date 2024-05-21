import numpy as np
from itertools import cycle

__all__ = ['ContinuumBasis', 'FourierBasis', 'PolynomialBasis']

class ContinuumBasis:
    
    def get_initial_guess(self, A, flux, ivar):
        ATCinv = A.T * ivar
        return np.linalg.solve(ATCinv @ A, ATCinv @ flux)
                

class FourierBasis(ContinuumBasis):    
        
    def __init__(self, P=7, L=None):
        self.P = P
        self.L = L
        return None
    
    @property
    def num_parameters(self):
        return self.P
    
    def get_or_set_L(self, λ):
        # We do this so that if we create a FourierBasis for fitting,
        # then later when we do the predictions we are using the same
        # length scale
        if self.L is None:
            self.L = 2 * np.ptp(λ)
        return self.L
    
    def design_matrix(self, λ):
        L = self.get_or_set_L(λ)            
        scale = (np.pi * λ) / L
        A = np.ones((λ.size, self.P), dtype=float)
        for j, f in zip(range(1, self.P), cycle((np.sin, np.cos))):
            A[:, j] = f(scale * (j + (j % 2)))        
        return A

            
class PolynomialBasis(ContinuumBasis):
    
    def __init__(self, deg=2):
        self.deg = deg
        return None
    
    @property
    def num_parameters(self):
        return self.deg + 1
    
    def design_matrix(self, λ):
        return np.vander(λ, self.deg + 1)
    
