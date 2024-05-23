
import numpy as np
import warnings
from scipy import stats, optimize as op
from types import MappingProxyType
from typing import Optional, Sequence, Tuple, Dict
from collections import OrderedDict

from Grok.models.linear_basis import LinearBasis, RadialBasis, PolynomialBasis, NeighbouringRadialBasis

    

# NOTES:
#   I think I would prefer to have a single SpectralLineModel and then provide some basis for nearby absorption
#   (which could be a radial basis or a NMF basis), but the problem is that when we are fitting with a radial
#   basis, we are fitting the line that we care about ALSO with the radial basis, whereas at the moment with
#   the NMF model we are fitting a profile. We could change this so that the NMF model ALSO uses a radial basis
#   to fit the EW of the line that we care about, OR we could change it so the RBF model isn't actually linear
#   any more and it uses a profile function to fit the line.

#   The NMF and RBF sets are just to represent the nearby absorption, so we should just separate that from the
#   main thing and REQUIRE that the central line absorption that we care about uses a RadialBasis


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
        
        if not isinstance(line_basis, LinearBasis):
            if line_basis is RadialBasis: # we know what to do here.                
                σ = 0.05 # default to something sensible.                
                offsets = np.linspace(-2, 2, 5)
                line_basis = RadialBasis(offsets * σ + λ, scales=σ)

            else:
                raise TypeError("`line_basis` must be an instance that is a sub-class of `LinearBasis`")
        
        if continuum_basis is not None and not isinstance(continuum_basis, LinearBasis):
            raise TypeError("`continuum_basis` must be an instance that is a sub-class of `LinearBasis`")
        
        if absorption_basis is not None and not isinstance(absorption_basis, LinearBasis):
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
        si, ei = spectrum.λ.searchsorted([self.λ - window, self.λ + window]) + [0, 1]
        λ, Y, Cinv = (spectrum.λ[si:ei], spectrum.flux[si:ei], np.diag(spectrum.ivar[si:ei]))
    
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
        fit = op.lsq_linear(ATCinv @ A + W, ATCinv @ Y, bounds.T) 
        
        if not fit.success:
            if raise_exception:
                raise RuntimeError(fit.message)
            else:
                warnings.warn(fit.message)
        
        # Get predictions from each base?
        base_predictions, si, y_pred = (OrderedDict(), 0, np.zeros_like(Y))
        for name, base in self.bases.items():
            sliced = slice(si, si + base.n_parameters)
            base_predictions[name] = A[:, sliced] @ fit.x[sliced]
            y_pred += base_predictions[name]
            si += base.n_parameters

        continuum = base_predictions.get("continuum")


        fig, ax = plt.subplots()
        ax.plot(λ, Y, c='k')
        ax.plot(λ, y_pred, c="tab:red")
        try:
            ax.axvspan(self.bases["absorption"].lower, self.bases["absorption"].upper, facecolor="#cccccc", zorder=-1)
        except:
            None
            
        colors = dict(line="tab:blue", absorption="tab:green", continuum="tab:orange")
        for name, pred in base_predictions.items():
            if name == "continuum":
                ax.plot(λ, pred, label=name, c=colors[name])            
            else:
                ax.plot(λ, pred + continuum, label=name, c=colors[name])
                
        for line in self.bases["line"].locations:
            ax.axvline(line, c="tab:blue", ls=":")
        ax.legend()
        ax.set_ylim(0.5, 1.1)

        
        
        
        

    
# NMF model needs:
# - wavelengths to compute at
# - spectral resolution to resample to
# - any shifts etc to apply to the stellar/telluric basis vectors
# then it can compute the design matrix for some wavelengths to 
    
    
if __name__ == "__main__":
    
    
    from Grok.specutils import Spectrum
    from Grok.models.continuum_basis import PolynomialBasis, FourierBasis
    
    spectrum = Spectrum.read("/Users/andycasey/software/smhr/smh/data/spectra/hd122563.fits")
    spectrum.ivar *= 1e4
    
    #model = RBFSpectralLineModel(5068.7655, continuum_basis=FourierBasis(3))
    #model.fit(spectrum, Lambda=1e8)

    import matplotlib.pyplot as plt
    
    wls = np.loadtxt("/Users/andycasey/Downloads/linelist_mm.txt", usecols=(0, 1, ))
    wls = wls[:, 0][wls[:, 1].astype(int) == 26]

    from Grok.models.linear_basis import NeighbouringGenericBasis 

    
    import h5py as h5
    with h5.File("Grok/bosz-highres-optical-basis-vectors.h5", "r") as fp:
        λ = fp["vacuum_wavelength"][:]
        stellar_basis_vectors = fp["stellar_basis_vectors"][:]
        meta = dict(fp["stellar_basis_vectors"].attrs)       
        try:
            telluric_basis_vectors = fp["telluric_basis_vectors"][:]
        except:
            telluric_basis_vectors = None        
    

    def vacuum_to_air(lambdas):
        """
        Convert vacuum wavelengths to air.

        The formula from Donald Morton (2000, ApJ. Suppl., 130, 403) is used for the 
        refraction index, which is also the IAU standard.

        As per https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion

        :param lambdas:
            The vacuum wavelengths.
        """
        s = 10**4 / lambdas
        n = 1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2)
        return (lambdas / n) 
    
    
    '''
    nmf_basis = NeighbouringGenericBasis(vacuum_to_air(λ), -stellar_basis_vectors.T, window=0)
    si, ei = spectrum.λ.searchsorted([wls[0], wls[-1]])    
    A = -nmf_basis.get_design_matrix(spectrum.λ[si:ei], 0)
    bounds = np.array([(0, np.inf)] * A.shape[1]).T    
    Y = 1 - spectrum.flux[si:ei]
    Cinv = np.eye(Y.size)
    lhs = A.T @ Cinv @ A + 1e-4 * np.eye(A.shape[1])
    rhs = A.T @ Cinv @ A @ Y
    result = op.lsq_linear(lhs, rhs, bounds)
    '''
    
    """
    fig, ax = plt.subplots()
    ax.plot(spectrum.λ[si:ei], spectrum.flux[si:ei], c='k')
    ax.plot(spectrum.λ[si:ei], 1 - A @ result.x, c="tab:red")
    """

    #nmf_basis = NeighbouringGenericBasis(vacuum_to_air(λ), -stellar_basis_vectors.T, X=result.x, X_tol=1e-2, window=0.5)

    from tqdm import tqdm
    for j, λ_line in enumerate(tqdm(wls)):

        #fig, axes = plt.subplots(2, 1, sharex=True, sharey=True)

        #model = RBFSpectralLineModel(λ_line, continuum_basis=PolynomialBasis(0))
        #model.fit(spectrum, Lambda=1e6, regularize_covariance=True, ax=axes[0])
        
        #model = RBFSpectralLineModel(λ_line, continuum_basis=PolynomialBasis(0))
        #model.fit(spectrum, ax=axes[1], Λ=1e4)
        '''
        spectrum.flux = 1-spectrum.flux
        model = SpectralLineModel(4714.393, line_basis=RadialBasis(4714.393, sign=1))
        model.fit(spectrum)
        '''
        
        absorption_basis = NeighbouringRadialBasis(λ_line - 0.5, λ_line + 0.5)
        
        model = LinearSpectralLineModel(λ_line, absorption_basis=absorption_basis)
        model.fit(spectrum)
        
        if j > 15:
            raise a