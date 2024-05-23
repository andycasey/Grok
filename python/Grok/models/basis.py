import h5py as h5
import numpy as np
from itertools import cycle
from scipy import stats
from typing import Union, Optional, Sequence, Dict

from Grok.specutils.lsf import instrument_lsf_kernel
from Grok.specutils.utils import apply_relativistic_velocity_shift
from Grok.utils import expand_path

__all__ = ['LinearBasis', 'FourierBasis', 'PolynomialBasis', 'RadialBasis']

class LinearBasis:
    
    """A base class for linear basis functions."""
    
    def get_bounds(self, bound=(-np.inf, np.inf)):
        return [bound] * self.n_parameters
    
    def _prepare_design_matrix(self, λ):
        self.design_matrix = self.get_design_matrix(λ)
    
    def get_initial_guess(self, λ, flux, ivar):
        try:
            A = self.design_matrix
        except AttributeError:            
            A = self.get_design_matrix(λ)
        ATCinv = A.T * ivar
        return np.linalg.solve(ATCinv @ A, ATCinv @ flux)
    
    def __call__(self, λ, *θ):
        try:
            return self.design_matrix @ θ
        except AttributeError:            
            return self.get_design_matrix(λ) @ θ



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
    
    """A linear basis for stellar spectra based on Non-negative Matrix Factorisation (NMF)."""
    
    def __init__(
        self,
        λ_rest_vacuum: Sequence[float],
        stellar_basis_vectors: Sequence[float],
        telluric_basis_vectors: Optional[Sequence[float]] = None,
        meta: Optional[Dict] = None
    ):
        """
        :param λ_rest_vacuum:
            An array of rest vacuum wavelengths (in Angstroms) of the basis vectors.
            The stellar basis vectors should be defined on this grid: at rest, in a vacuum.
            Ironically, for consistency we also require the telluric basis vectors to be defined 
            on this grid: at rest, in a vacuum.
        
        :param stellar_basis_vectors:
            A (C, P)-shape array of stellar basis vectors with non-negative entries. Here, C is
            the number of basis vectors and P is the same as the size of `λ_rest_vacuum`.
            
            These basis vectors should be defined to be in the rest frame, in a vacuum.
                    
        :param telluric_basis_vectors: [optional]
            A (B, P)-shape array of telluric basis vectors with non-negative entries. Here, P is
            the same as the size of `λ_rest_vacuum`.
            
            Ironically, even though tellurics are defined in the observed frame (usually in air),
            these telluric basis vectors must be defined in the rest frame, in a vacuum.
            
        :param meta: [optional]
            A metadata dictionary.            
        """
        self.λ_rest_vacuum = λ_rest_vacuum
        self.stellar_basis_vectors = stellar_basis_vectors
        self.telluric_basis_vectors = telluric_basis_vectors
        self.meta = meta or dict()
        return None        
    
    
    @property
    def n_parameters(self):
        """The number of parameters in the model."""
        n = self.stellar_basis_vectors.shape[0]
        if self.telluric_basis_vectors is not None:
            n += self.telluric_basis_vectors.shape[0]
        return n
    
    
    def get_bounds(self, bound=(0, np.inf)):
        """Return the bounds on the model parameters."""
        return [bound] * self.n_parameters
    
    def get_design_matrix(
        self,
        λ_observed_vacuum: Sequence[float],
        *args,
        stellar_v_rel: Optional[float] = None,
        telluric_v_rel: Optional[float] = None,
        Ro: Optional[float] = None,
        Ri: Optional[float] = None,
        **kwargs
    ):
        """
        Get the design matrix for the given observed wavelengths.
        
        :param λ_observed_vacuum:
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
        
        use_convolution = (Ro is not None)
            
        # shift this model to the stellar frame, then interpolate to the observed points
        λ_stellar_vacuum = apply_relativistic_velocity_shift(self.λ_rest_vacuum, stellar_v_rel or 0)
        
        if use_convolution:
            Ri = Ri or self.meta.get("Ri", np.inf)
            if Ro > Ri:
                raise ValueError("The output spectral resolution must be lower than the input spectral resolution.")
            
            # The input wavelengths are vacuum stellar rest frame, and output are observed vacuum.
            K_stellar = instrument_lsf_kernel(λ_stellar_vacuum, λ_observed_vacuum, Ro, Ri)
            stellar_basis_vectors = self.stellar_basis_vectors @ K_stellar

            if self.telluric_basis_vectors is not None:
                if stellar_v_rel == telluric_v_rel:
                    # Special (improbable) case where the stellar and telluric velocities are the same.
                    K_telluric = K_stellar
                else:
                    K_telluric = instrument_lsf_kernel(
                        apply_relativistic_velocity_shift(self.λ_rest_vacuum, telluric_v_rel or 0),
                        λ_observed_vacuum,
                        Ro,
                        Ri
                    )
                telluric_basis_vectors = self.telluric_basis_vectors @ K_telluric
            else:
                telluric_basis_vectors = np.empty((0, λ_observed_vacuum.size))
            
        else:        
            # Interpolation only.            
            stellar_basis_vectors = _interpolate_basis_vector(λ_observed_vacuum, λ_stellar_vacuum, self.stellar_basis_vectors)            
            if self.telluric_basis_vectors is not None:                
                telluric_basis_vectors = _interpolate_basis_vector(
                    λ_observed_vacuum, 
                    apply_relativistic_velocity_shift(self.λ_rest_vacuum, telluric_v_rel or 0),
                    self.telluric_basis_vectors
                )
            else:
                telluric_basis_vectors = np.empty((0, λ_observed_vacuum.size))
    
        A = np.vstack([stellar_basis_vectors, telluric_basis_vectors]).T
        return A


    @classmethod
    def from_path(cls, path):
        """
        Load a model from disk.
        
        :param path:
            The path to the saved model.
        """
        with h5.File(expand_path(path), "r") as fp:
            λ_vacuum = fp["vacuum_wavelength"][:]
            stellar_basis_vectors = fp["stellar_basis_vectors"][:]
            meta = dict(fp["stellar_basis_vectors"].attrs)       
            try:
                telluric_basis_vectors = fp["telluric_basis_vectors"][:]
            except:
                telluric_basis_vectors = None
                
        return cls(
            λ_vacuum,
            stellar_basis_vectors,
            telluric_basis_vectors=telluric_basis_vectors,
            meta=meta,            
        )
    
    
    
def _interpolate_basis_vector(λ_output, λ_input, basis_vectors, left=0, right=0):
    bv = np.zeros((basis_vectors.shape[0], λ_output.size))
    for c, basis_vector in enumerate(basis_vectors):
        bv[c] = np.interp(λ_output, λ_input, basis_vector, left=left, right=right)
    return bv    