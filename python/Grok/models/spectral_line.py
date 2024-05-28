
import numpy as np
import warnings
from scipy import optimize as op
from types import MappingProxyType
from typing import Optional, Dict, Sequence, Union

from Grok.models.basis import LinearBasis, RadialBasis, PolynomialBasis
from Grok.models.profile import Profile, GaussianProfile
from Grok.specutils.utils import ivar_to_sigma


class SpectralLineModel:
    
    """Model an absorption line in a spectrum."""
    
    def __init__(
        self,
        λ,
        profile: Optional[Profile] = GaussianProfile,
        bases: Optional[Union[LinearBasis, Dict[str, LinearBasis]]] = None,
        mask_regions: Optional[Sequence[Sequence[float]]] = None
    ):
        """
        Model an absorption line in a spectrum.
        
        :param λ:
            The central wavelength of the absorption line.
            
        :param profile: [optional]
            The profile to model the absorption line.
        
        :param bases: [optional]
            Some optional bases to model the continuum or nearby absorption features. Bases enter multiplicatively.
            
        :param mask_regions: [optional]
            A list of regions to exclude from the fit. This should be a two-length tuple for each region.
        """
        # TODO: Consider whether to include an absorption_basis option here
        self.λ = λ
        self.profile = profile(λ) if isinstance(profile, type) else profile
        if bases is None:
            self.bases = MappingProxyType({})
        elif isinstance(bases, LinearBasis):
            self.bases = MappingProxyType({None: bases})
        else:
            self.bases = MappingProxyType(bases)            
        self.mask_regions = mask_regions
        return None

    @property
    def n_parameters(self):
        """The total number of model parameters."""
        return self.profile.n_parameters + sum([b.n_parameters for name, b in self.bases.items()])
    
    
    def get_mask(self, λ):
        """
        Get a mask for the regions to exclude from the fit.
        
        :param λ:
            The wavelengths to fit.
        """
        mask = np.zeros(λ.size, dtype=bool)
        if self.mask_regions is not None:
            for region in self.mask_regions:
                mask |= (region[1] >= λ) & (λ >= region[0])            
        return mask
            
    
    def prepare_fit(self, spectrum, window: Optional[float] = 2):
        """
        Given a spectrum and a window to fit around, prepare the data and design matrices.
        
        :param spectrum:
            The spectrum to fit.
            
        :param window:
            The window to fit around the central wavelength. The full fitting domain is twice this window.
        """
        # Get the data.        
        si, ei = spectrum.λ_rest_vacuum.searchsorted([self.λ - window, self.λ + window]) + [0, 1]
        λ, flux, ivar = (spectrum.λ_rest_vacuum[si:ei], spectrum.flux[si:ei], spectrum.ivar[si:ei])
    
        use = ~self.get_mask(λ)
            
        # Build design matrices and get bounds.
        A, bounds = ([], self.profile.get_bounds())
        p0 = list(self.profile.get_initial_guess(λ[use], flux[use], ivar[use]))
        for name, base in self.bases.items():
            A.append(base.get_design_matrix(λ[use], self.λ))
            p0.extend(base.get_initial_guess(λ[use], flux[use]))
            bounds.extend(base.get_bounds())            
        return (λ, flux, ivar, use, A, bounds, p0)       


    def get_bounds(self, *args, **kwargs):
        bounds = self.profile.get_bounds()
        for base in self.bases.values():
            bounds.extend(base.get_bounds())
        return bounds


    def fit(
        self,
        spectrum,
        p0: Optional[Sequence[float]] = None,
        window: Optional[float] = 2,
        **kwargs
    ):    
        """
        Fit the absorption line model to a spectrum.
        
        :param spectrum:
            The spectrum to fit.
            
        :param p0: [optional]
            The initial guess for the model parameters.
        
        :param window: [optional]
            The window to fit around the central wavelength. The full fitting domain is twice this window.
        """
        λ, flux, ivar, use, As, bounds, x0 = self.prepare_fit(spectrum, window)
        
        n = self.profile.n_parameters
            
        def f(λ, *θ):
            y = self.profile(λ, *θ[:n])
            si = n
            for A, base in zip(As, self.bases.values()):
                y *= base(λ, θ[si:si+A.shape[1]], A=A)
                si += A.shape[1]
            return y
        
        if p0 is None:
            p0 = x0
                
        self.θ, self.Σ = op.curve_fit(
            f,
            λ[use],
            flux[use],
            p0=p0,
            sigma=ivar_to_sigma(ivar[use]),
            bounds=np.array(bounds).T
        )
        si = 0
        for A, base in zip(As, self.bases.values()):
            n = A.shape[1]
            base.θ = self.θ[si:si+n]
            base.Σ = self.Σ[si:si+n, si:si+n]
            si += n        
        return self
    
    
    @property
    def equivalent_width(self):
        try:
            self.θ
        except AttributeError:
            return None
        else:
            return self.profile.get_equivalent_width(self.θ)
        
    @property
    def equivalent_width_uncertainty(self):
        try:
            args = (self.θ, self.Σ)
        except AttributeError:
            return None
        else:
            return self.profile.get_equivalent_width_uncertainty(*args)
    
    
    def __call__(self, λ, *θ, full_output=False):
        θ = self.θ if len(θ) == 0 else θ
        n = self.profile.n_parameters
        y = self.profile(λ, *θ[:n])

        y_bases = dict(profile=np.copy(y))
        si = n
        for name, base in self.bases.items():
            A = base.get_design_matrix(λ, self.λ)
            y_bases[name] = base(λ, θ[si:si + base.n_parameters], A=A)
            y *= y_bases[name]
            si += base.n_parameters
        
        if full_output:
            return (y, y_bases)
        return y
        
    
    def plot(self, spectrum, ax=None, **kwargs):        
        λ, flux, *_ = self.prepare_fit(spectrum, **kwargs)
        
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        
        y, y_bases = self(λ, full_output=True)
        ax.plot(λ, flux, c='k', label="Data")
        ax.plot(λ, y, c="tab:red", label="Model")
        for name, pred in y_bases.items():
            ax.plot(λ, pred, label=name)
            
        ax.axvline(self.λ, c="#666666", ls=":", zorder=-1, lw=0.5)
        if self.mask_regions is not None:
            for s, e in self.mask_regions:
                ax.axvspan(s, e, 0, 1.2, facecolor="#CCCCCC", zorder=-1, alpha=0.5)
                
        try:
            v = (
                self.profile.λ - self.profile.λ_tolerance,
                self.profile.λ + self.profile.λ_tolerance
            )
        except:
            None
        else:
            for _ in v:
                ax.axvline(_, c="#666666", ls="-.", zorder=-1, lw=0.5)
        
        ax.set_xlim(λ[[0, -1]])
        ax.legend()
        ax.set_xlabel(r"Wavelength $\lambda$ [vacuum; rest]")
        ax.set_ylabel("Flux")
        
        return fig



class LinearSpectralLineModel:
    
    """Model an absorption line in a spectrum as a sum of linear components."""
    
    def __init__(
        self,
        λ,
        line_basis: Optional[LinearBasis] = RadialBasis,
        continuum_basis: Optional[LinearBasis] = PolynomialBasis(0),
        absorption_basis: Optional[LinearBasis] = None,
    ):
        """
        Model an absorption line in a spectrum.
        
        :param λ: 
            The central wavelength of the absorption line.
        
        :param line_basis: [optional]
            The basis to model the absorption line. 
            
        :param continuum_basis: [optional]
            The basis to model the continuum. 
        
        :param absorption_basis: [optional]
            The basis to model nearby absorption features in the spectrum. 
        """
        
        # TODO: consider whether to include a `mask_regions` here.
        
        if not issubclass(line_basis, LinearBasis):
            raise TypeError("`line_basis` must be an instance that is a sub-class of `LinearBasis`")
        else:
            if line_basis is RadialBasis: # we know what to do here.                
                σ = 0.05 # default to something sensible MAGIC 
                offsets = np.linspace(-2, 2, 5)
                line_basis = RadialBasis(offsets * σ + λ, scales=σ)
            # if it's not instantiated, we don't know what to do with it
        
        if continuum_basis is not None and isinstance(continuum_basis, type) and not issubclass(continuum_basis, LinearBasis):
            raise TypeError("`continuum_basis` must be an instance that is a sub-class of `LinearBasis`")
        
        if absorption_basis is not None and isinstance(continuum_basis, type) and not issubclass(absorption_basis, LinearBasis):
            raise TypeError("`absorption_basis` must be an instance that is a sub-class of `LinearBasis`")
        
        all_bases = [
            ("line", line_basis),
            ("absorption", absorption_basis),
            ("continuum", continuum_basis)
        ]                
        self.λ = λ
        self.bases = MappingProxyType({ label: base for label, base in all_bases if base is not None })
        return None


    @property
    def n_parameters(self):
        """The total number of model parameters."""
        return sum(v.n_parameters for v in self.bases.values())


    def prepare_fit(self, spectrum, window: Optional[float] = 2):
        """
        Given a spectrum and a window to fit around, prepare the data and design matrices.
        
        :param spectrum:
            The spectrum to fit.
        
        :param window:
            The window to fit around the central wavelength. The full fitting domain is twice this window.
        """
        
        # Get the data.        
        si, ei = spectrum.λ_rest_vacuum.searchsorted([self.λ - window, self.λ + window]) + [0, 1]
        λ, Y, Cinv = (spectrum.λ_rest_vacuum[si:ei], spectrum.flux[si:ei], np.diag(spectrum.ivar[si:ei]))
    
        # Build design matrices and get bounds.
        args, bounds = ([], [])
        for base in self.bases.values():
            args.append(base.get_design_matrix(λ, self.λ))
            bounds.extend(base.get_bounds())
        
        A = np.hstack(args)
        bounds = np.atleast_2d(bounds)

        return (λ, Y, Cinv, A, bounds)
    
    
    def get_feature_weight_matrix(self, Λ):
        """
        Given a dictionary of regularization strengths, return a diagonal matrix of feature weights.
        
        :param Λ:
            A dictionary of regularization strengths for each base, with the base name as the dictionary key,
            and the regularization strength(s) as the value. The regularization strength can be a float, or
            an `n`-length array, where `n` is the number of parameters in that base.
        """
        Λ = Λ.copy()
        W = np.hstack([Λ.pop(k, 0) * np.ones(v.n_parameters) for k, v in self.bases.items()])
        unknown_bases = set(Λ).difference({"line", "absorption", "continuum"})
        if len(unknown_bases) > 0:
            warnings.warn(f"Ignoring regularization strengths for unknown base(s): {', '.join(unknown_bases)}")
        return np.diag(W)
    
    
    def __call__(self, λ, *θ):        
        θ = self.θ if len(θ) == 0 else θ        
        A = np.hstack([base.get_design_matrix(λ, self.λ) for base in self.bases.values()])
        return A @ θ
            
        
    def fit(
        self,
        spectrum,
        window: Optional[float] = 2,
        Λ: Optional[Dict[str, float]] = MappingProxyType(dict(line=0, absorption=0, continuum=0)),
        raise_exception: Optional[bool] = False,
        **kwargs
    ):
        """
        Fit the absorption line model to a spectrum.
        
        :param spectrum:
            The spectrum to fit.
            
        :param window: [optional]
            The window to fit around the central wavelength. The full fitting domain is twice this window.
        
        :param Λ: [optional]
            A dictionary of regularization strengths for each base, with the base name as the dictionary key,
            and the regularization strength as the value. The regularization strength can be a float, or
            an `n`-length array, where `n` is the number of parameters in that base.
        
        :param raise_exception: [optional]
            Raise an exception if the optimization fails.
        """
        
        λ, Y, Cinv, A, bounds = self.prepare_fit(spectrum, window)
        
        W = self.get_feature_weight_matrix(Λ)
        
        ATCinv = A.T @ Cinv        
        self.fit = op.lsq_linear(ATCinv @ A + W, ATCinv @ Y, bounds.T) 
        
        if not self.fit.success:
            if raise_exception:
                raise RuntimeError(self.fit.message)
            else:
                warnings.warn(self.fit.message)
        
        self.θ = self.fit.x
        return None
    
    
    @property
    def equivalent_width(self):
        raise NotImplementedError("needs a thinko")
    
