import h5py as h5
import numpy as np
from itertools import cycle
from scipy import stats
from typing import Union, Optional, Sequence, Dict
from types import MappingProxyType
from functools import cached_property
import warnings

from Grok.spectrum.lsf import instrument_lsf_kernel
from Grok.spectrum.utils import apply_relativistic_velocity_shift
from Grok.utils import expand_path

__all__ = ['LinearBasis', 'FourierBasis', 'PolynomialBasis', 'RadialBasis']

class LinearBasis:
    
    """A base class for linear basis functions."""
    
    def get_bounds(self, bound=(-np.inf, np.inf)):
        return [bound] * self.n_parameters
        
    def get_initial_guess(self, λ, flux, ivar=None, **kwargs):
        A = kwargs.get("A", None)
        if A is None:
            A = self.get_design_matrix(λ)
        
        if ivar is None:
            ivar = 1
        ATCinv = A.T * ivar
        return np.linalg.solve(ATCinv @ A, ATCinv @ flux)
    
    def __call__(self, λ, θ, **kwargs):
        A = kwargs.get("A", None)
        if A is None:
            A = self.get_design_matrix(λ)
        
        return np.dot(A, θ)
    
    @property
    def n_parameters(self):
        raise ValueError("Subclasses must implement this property.")


class FourierBasis(LinearBasis):
    
    """A set of sine-and-cosine basis functions."""
        
    def __init__(self, P=7, L=None):
        self.P = P
        self.L = L
        return None
    
    @property
    def n_parameters(self):
        return self.P
    
    def get_or_set_L(self, λ):
        # We do this so that if we create a FourierBasis for fitting,
        # then later when we do the predictions we are using the same
        # length scale
        if self.L is None:
            self.L = 2 * np.ptp(λ)
        return self.L
    
    def get_design_matrix(self, λ, *args, **kwargs):
        L = self.get_or_set_L(λ)            
        scale = (np.pi * λ) / L
        A = np.ones((λ.size, self.P), dtype=float)
        for j, f in zip(range(1, self.P), cycle((np.sin, np.cos))):
            A[:, j] = f(scale * (j + (j % 2)))        
        return A
    

            
class PolynomialBasis(LinearBasis):
    
    """A polynomial basis set."""
    
    def __init__(self, deg=2):
        self.deg = deg
        return None
    
    @property
    def n_parameters(self):
        return self.deg + 1
    
    def get_design_matrix(self, λ, *args, **kwargs):
        return np.vander(λ, self.deg + 1)
    


class RadialBasis(LinearBasis):
    
    def __init__(self, locations, scales=0.05, sign=-1):
        self.locations = np.atleast_1d(locations).flatten()
        if isinstance(scales, (int, float)):
            self.scales = np.ones_like(self.locations) * scales
        else:
            self.scales = scales
        self.sign = np.sign(sign)
        return None
    
    @property
    def n_parameters(self):
        return self.locations.size
                
    def get_design_matrix(self, λ, *args, **kwargs):
        A = np.zeros((λ.size, self.n_parameters), dtype=float)
        for i, (loc, scale) in enumerate(zip(self.locations, self.scales)):
            A[:, i] = np.sign(self.sign) * stats.norm.pdf(λ, loc=loc, scale=scale)        
        return A

    def get_bounds(self, bound=(0, np.inf)):
        return [bound] * self.n_parameters


class NeighbouringGenericBasis(LinearBasis):
    
    def __init__(self, λ, components, window, X=None, X_tol=1e-6):
        self.λ = λ
        self.components = np.atleast_2d(components)
        assert self.λ.size == self.components.shape[0]
        if X is None:
            self._bounds = [(0, np.inf)] * self.n_parameters
        else:
            self._bounds = np.clip(
                np.array([X - X_tol, X + X_tol]),
                0, 
                np.inf
            ).T
        self.window = window
        return None
    
    @property
    def n_parameters(self):
        return self.components.shape[1]
    
    def get_bounds(self, **kwargs):
        return self._bounds
        
    def get_design_matrix(self, λ, λ_line, *args, **kwargs):
        # interpolate
        A = np.zeros((λ.size, self.n_parameters), dtype=float)
        for i in range(self.n_parameters):
            A[:, i] = np.interp(λ, self.λ, self.components[:, i], left=0, right=0)
        
        # central mask
        A[(λ >= (λ_line - self.window)) * (λ <= (λ_line + self.window))] = 0
        
        return A
    

class NeighbouringRadialBasis(LinearBasis):
    
    """A set of radial basis functions that neighbour a masked region (e.g., a line to measure)."""
    
    def __init__(self,
        lower: float,
        upper: float,
        locations: Optional[Sequence[float]] = None,
        scales: Optional[Union[float, Sequence[float]]] = 0.05,
        n: Optional[int] = 100,
        sign: Optional[float] = -1
    ):
        """
        Define a basis that includes radial basis functions across the domain to be 
        fit, except for between the `lower` and `upper` positions. 
        
        If no specific `locations` are given, then the bases will be placed equispaced 
        throughout the domain to be fit.
        
        :param lower:
            The lower mask position. No radial basis functions will be placed between
            the `lower` and `upper` positions, but if the `scales` are large enough,
            there might be support in this region from nearby bases.
        
        :param upper:
            The upper mask position. No radial basis functions will be placed between
            the `lower` and `upper` positions, but if the `scales` are large enough,
            there might be support in this region from nearby bases.
        
        :param locations: [optional]
            The locations where to place radial basis functions. If `None` is given
            
        :param scales: [optional]
            The scales of the radial basis functions.
        
        :param n: [optional]
            The number of radial basis functions to use. If `locations` is None,
            then this will be the number of radial basis functions placed throughout
            the domain.
        """
        self.lower = lower
        self.upper = upper
        self.locations = locations
        self.scales = scales
        self.n = n
        self.sign = np.sign(sign)
        
        if self.locations is not None:
            self.locations = np.array(self.locations).flatten()            
            if np.any((upper >= self.locations) * (self.locations >= lower)):
                raise ValueError("Some locations are within the mask region.")
                
        if (
            self.locations is not None 
        and np.atleast_1d(self.locations).size != np.atleast_1d(self.scales).size
        ):
            raise ValueError("Number of scales does not match number of locations.")

        if self.n % 2 > 0:
            raise ValueError("`n` must be an even number.")
        return None
    
    
    def get_design_matrix(self, λ, *args, **kwargs):
        A = np.zeros((λ.size, self.n), dtype=float)        
        if self.locations is None:
            # TODO: check these locations match my intuition
            self.locations = np.hstack([        
                np.linspace(λ[0], self.lower, self.n // 2 + 1)[:-1],
                np.linspace(self.upper, λ[-1], self.n // 2 + 1)[1:]
            ])
                
        for i, (μ, σ) in enumerate(zip(self.locations, cycle(np.atleast_1d(self.scales)))):            
            A[:, i] = self.sign * stats.norm.pdf(λ, μ, σ)
        return A
            
    @property
    def n_parameters(self):
        return self.locations.size
    
    def get_bounds(self, bound=(0, np.inf)):
        return [bound] * self.n_parameters
    
    
    
class NMFSpectralBasis(LinearBasis):
    
    """A linear basis for spectra based on Non-negative Matrix Factorisation (NMF)."""
    
    def __init__(
        self,
        λ_rest_vacuum: Sequence[float],
        basis_vectors: Sequence[float],
        mask_regions: Optional[Sequence[float]] = None,
        v_rel: Optional[float] = None,
        Ro: Optional[float] = None,
        Ri: Optional[float] = None,
        **kwargs,
    ):
        """
        :param λ_rest_vacuum:
            An array of rest vacuum wavelengths (in Angstroms) of the basis vectors.
            The stellar basis vectors should be defined on this grid: at rest, in a vacuum.
            Ironically, for consistency we also require the telluric basis vectors to be defined 
            on this grid: at rest, in a vacuum.
        
        :param basis_vectors:
            A (C, )-shape array of basis vectors with non-negative entries. Here, C is
            the number of basis vectors and P is the same as the size of `λ_rest_vacuum`.
            
            These basis vectors should be defined to be in the rest frame, in a vacuum.
        
        :param mask_regions: [optional]
            A list of tuples, each defining a mask region. Each tuple should have two elements:
            the lower and upper bounds of the mask region.
            
        :param v_rel: [optional]
            The relative velocity of the observer with respect to the rest frame of the source.
        
        :param Ro: [optional]
            The spectral resolution of the observations. If `None` is given, then no convolution
            will take place; the basis vectors will be interpolated.
        
        :param Ri: [optional]
            The spectral resolution of the basis vectors. If `None` is given then it defaults to
            infinity.        
        """
        self.λ_rest_vacuum = λ_rest_vacuum
        self.basis_vectors = basis_vectors        
        self.mask_regions = mask_regions
        if self.basis_vectors.shape[0] < self.basis_vectors.shape[1]:
            raise ValueError("I don't believe you have more bases than pixels.")
        self.taper_region = None
        self.v_rel = v_rel
        self.Ro = Ro
        self.Ri = Ri
        if Ro is not None and Ri is not None and Ro > Ri:
            raise ValueError("The output spectral resolution must be lower than the input spectral resolution.")        
        return None   


    def __copy__(self):
        """Return a copy of this basis."""
        k = self.__class__(
            self.λ_rest_vacuum,
            self.basis_vectors,
            self.mask_regions,
            self.v_rel,
            self.Ro,
            self.Ri
        )
        try:
            k.θ = self.θ
        except AttributeError:
            None
        return k
    

    def set_taper_region(self, lower, upper, scale=0.05):
        self.taper_region = (lower, upper, scale)
        return self


    @property
    def n_parameters(self):
        """The number of parameters in the model."""
        return self.basis_vectors.shape[1]


    def get_bounds(self, bound=(0, np.inf)):
        """Return the bounds on the model parameters."""
        return [bound] * self.n_parameters
    
        
    def get_initial_guess(self, *args, **kwargs):
        try:
            return self.θ
        except AttributeError:
            return 1e-5 * np.ones(self.n_parameters)
    
    
    def get_design_matrix(self, λ: Sequence[float], *args, **kwargs):
        """
        Get the design matrix for the given observed wavelengths.
        
        :param λ:
            The observed vacuum wavelengths to compute this design matrix for.
        
        :param stellar_v_rel: [optional]
            The relative velocity of the observer with respect to the rest frame of the
            source. If `None` is given, then no radial velocity shift will be applied.
                
        :param Ro: [optional]
            The spectral resolution of the observations. If `None` is given, then no
            convolution will take place; the basis vectors will be interpolated.
        
        :param Ri: [optional]
            The spectral resolution of the basis vectors. If `None` is given then
            it defaults to the input spectral resolution stored in the metadata of
            the basis vector file, or infinity if the input spectral resolution is
            not stored in that file.
        """
                    
        # shift this model to the stellar frame, then interpolate to the observed points
        λ_shift_vacuum = apply_relativistic_velocity_shift(self.λ_rest_vacuum, self.v_rel or 0)
        
        if self.Ro is not None: # Do convolution            
            # The input wavelengths are vacuum stellar rest frame, and output are observed vacuum.
            kernel = kwargs.pop("kernel", None)
            if kernel is None:
                kernel = instrument_lsf_kernel(λ_shift_vacuum, λ, self.Ro, self.Ri or np.inf)
            A = (self.basis_vectors.T @ kernel).T
        else:        
            # Interpolation only.            
            A = _interpolate_basis_vector(λ, λ_shift_vacuum, self.basis_vectors)
                
        if self.mask_regions is not None:
            # We don't model things in the mask region.
            for lower, upper in self.mask_regions:
                A[(λ >= lower) * (λ <= upper)] = 0
        
        if self.taper_region is not None:
            lower, upper, scale = self.taper_region
            
            si, ei = λ.searchsorted([lower, upper]) + [0, 1]
            xi = λ[si:ei]
            yi = np.ones_like(xi)
            lhs = stats.norm.pdf(xi[:xi.size // 2], lower, scale)
            yi[:xi.size // 2] = lhs / lhs[0]
            rhs = stats.norm.pdf(xi[xi.size // 2:], upper, scale)
            yi[xi.size // 2:] = rhs / rhs[-1]
            A[si:ei] *= yi[:, None]
        return A


    def __call__(self, λ, θ, **kwargs):
        A = kwargs.get("A", None)
        if A is None:
            A = self.get_design_matrix(λ)        
        return 1 - np.dot(A, θ)
    
    
    
def _interpolate_basis_vector(λ_output, λ_input, basis_vectors, left=0, right=0):
    bv = np.zeros((λ_output.size, basis_vectors.shape[1]))
    for c, basis_vector in enumerate(basis_vectors.T):
        bv[:, c] = np.interp(λ_output, λ_input, basis_vector, left=left, right=right)
    return bv    
