
import numpy as np
import warnings
from scipy import stats, optimize as op
from typing import Optional, Sequence, Tuple



from Grok.models.continuum_basis import ContinuumBasis, RadialBasis, PolynomialBasis


class BaseSpectralLineModel:
    pass
    


class RBFSpectralLineModel(BaseSpectralLineModel):
    
    """
    Fit an absorption line and the surrounding absorption features using a radial basis function model.
    """
    
    def __init__(
        self,
        λ,
        rbf_centroids: Optional[Sequence[float]] = None,
        n_rbf: int = 101,
        sigma: Optional[float] = 0.05,
        window: Optional[float] = 2,
        continuum_basis: Optional[ContinuumBasis] = None,
        mask_regions: Optional[Sequence[Tuple[float, float]]] = None,
    ):
        self.λ = λ
        self.sigma = sigma
        self.window = window
        self.continuum_basis = continuum_basis
        self.mask_regions = mask_regions
        
        if rbf_centroids is None:        
            self.n_rbf = n_rbf
            self.rbf_centroids = None
            if (self.n_rbf % 2) == 0:
                raise ValueError("Number of RBF basis must be odd, or supply your own `rbf_centroids`.")                    
        elif rbf_centroids is not None and n_rbf is not None:
            raise ValueError("Supply either `rbf_centroids` or `n_rbf`, not both.")
        else:
            self.rbf_centroids = np.array(rbf_centroids)
            self.n_rbf = len(self.rbf_centroids)
            
        return None
    
    @property
    def n_parameters(self):
        n = self.n_rbf 
        try:
            n += self.continuum_basis.n_parameters
        except AttributeError:
            None
        finally:
            return n

    
    def _prepare_fit(self, spectrum):
        # Get the data.        
        si, ei = spectrum.λ.searchsorted([self.λ - self.window, self.λ + self.window]) + [0, 1]
        λ, Y, Cinv = (spectrum.λ[si:ei], spectrum.flux[si:ei], np.diag(spectrum.ivar[si:ei]))
    
        if self.rbf_centroids is None:            
            # Build design matrices.    
            rbf_centroids = np.hstack([
                np.linspace(λ[0], self.λ, self.n_rbf // 2 + 1)[:-1],
                np.linspace(self.λ, λ[-1], self.n_rbf // 2 + 1)            
            ])
        else:
            rbf_centroids = self.rbf_centroids
            
        rb = RadialBasis(rbf_centroids, self.sigma)
        A_rb = rb.get_design_matrix(λ)
        bounds = [(0, np.inf)] * rb.n_parameters
        if self.continuum_basis is not None:
            A_c = self.continuum_basis.get_design_matrix(λ)
            bounds.extend(self.continuum_basis.get_bounds())
        else:
            A_c = np.empty((λ.size, 0))
        
        A = np.hstack([A_rb, A_c])
        return (λ, Y, Cinv, rb, A_rb, A_c, A, np.atleast_2d(bounds).T)
                
                
    def _fit_least_squares(self, A, flux, ivar, Lambda, full_output=False):
        # Fit the model.
        ATCinv = A.T * ivar
        ATCinvA = ATCinv @ A
        ATCinvY = ATCinv @ flux
        R = Lambda * np.ones(A.shape[1])
        R[self.n_rbf:] = 0 # do not regularize the continuum coefficients
        central = 3
        R[self.n_rbf // 2 - central:self.n_rbf // 2 + central] = 0 # do not regularize the central few bases 
        R = np.diag(R)

        θ = np.linalg.solve(ATCinvA + R, ATCinvY)
        
        if full_output:
            return (θ, ATCinv, ATCinvA, ATCinvY)
        return θ
            
    
    def get_feature_weight_matrix(self, Λ, ignore_central_rbf=3, **kwargs):
        W = Λ * np.ones(self.n_parameters)
        W[self.n_rbf:] = 0 # don't regularize continuum coefficients
        W[self.n_rbf // 2 - ignore_central_rbf:self.n_rbf // 2 + ignore_central_rbf] = 0 # don't regularize central RBF kernels
        return np.diag(W)
    
        
    def fit(
        self,
        spectrum,
        Λ: Optional[float] = 1e6,
        ax=None,
        raise_exception=False,
        **kwargs
    ):
        λ, Y, Cinv, locations, A_rb, A_c, A, bounds = self._prepare_fit(spectrum)
        
        W = self.get_feature_weight_matrix(Λ, **kwargs)
        
        ATCinv = A.T @ Cinv        
        fit = op.lsq_linear(ATCinv @ A + W, ATCinv @ Y, bounds) 
        
        if not fit.success:
            if raise_exception:
                raise RuntimeError(fit.message)
            else:
                warnings.warn(fit.message)
        
        # Get the continuum and line components.
        pred_absorption = 
        pred_continuum = A_c @ fit.x[self.n_rbf:]
        
        
        central = 3

        if ax is None:
            fig, ax = plt.subplots()
            
        
        ax.plot(λ, flux, c='k')
        ax.plot(λ, A @ θ, c="tab:blue")
        
        central_mask = np.zeros(θ.size, dtype=bool)
        central_mask[self.n_rbf // 2 - central:self.n_rbf // 2 + central] = True
        central_mask[self.n_rbf:] = True
        ax.plot(λ, A[:, central_mask] @ θ[central_mask], c="tab:green")
        
        ax.set_ylim(0.5, 1.1)
        
        # find contiguous places
        is_non_zero = (θ[:self.n_rbf] > (1e-3)) # anything more than 1 mA
        edges = 1 + np.where(np.diff(is_non_zero))[0]
        indices = np.sort(np.hstack([
            0, np.repeat(edges, 2), self.n_rbf 
        ]).reshape((-1, 2)))
        
        non_zero_indices = indices[0 if is_non_zero[0] else 1::2]
        for si, ei in non_zero_indices:
            assert is_non_zero[si:ei].all()
            assert si != ei
            v = np.sum(locations[si:ei] * θ[si:ei] / np.sum(θ[si:ei]))
            ew = 1000 * np.sum(θ[si:ei])
            if ew > 5:
                print(f"{si} {ei} {v:.3f} {ew:.3f}")
                ax.axvline(v, c="#666666", ls=":", zorder=-1)
                ax.axvspan(locations[si], locations[ei-1], color="#cccccc", alpha=0.5, zorder=-1)

        #for loc in locations:
        #    ax.axvline(loc, c="tab:red", lw=0.5, zorder=-1)
        
        ax.axvline(λ_line, c="tab:blue")
        print(result)
        print(θ)
        return None
        
    
    def old_code_for_fit(
        self,
        spectrum,
        op_kwds: Optional[dict] = None,
        Lambda: Optional[float] = 1e3,
        regularize_covariance: Optional[bool] = False,
        ax=None,
        **kwargs
    ):
        λ, flux, ivar, locations, A_rb, A_c, A, bounds = self._prepare_fit(spectrum)

        kwds = dict(max_iter=10_000)
        kwds.update(op_kwds or {})
        
        # Let's get a bounded unregularized estimate that does not account for uncertainties.
        # This seemed like a good idea but was sometimes bad in practice. 
        # going to use he least squares result instead, even though it is unbounded and a different regularization scheme (L1 vs L2)
        #p_init = op.lsq_linear(A, flux, bounds, **kwds)        
        # Also, it turns out lsq_linear is very slow...

        x0 = self._fit_least_squares(A, flux, ivar, Lambda)

        # Create a weighting matrix for regularization
        central = 3 # don't regularize the central few bases (either side)
        W = np.sqrt(Lambda) * np.ones_like(x0)
        W[self.n_rbf:] = 0 # don't regularize continuum coefficients
        W[self.n_rbf // 2 - central:self.n_rbf // 2 + central] = 0 # don't regularize central RBF kernels
        
        corrcoef = np.eye(W.size)
        if regularize_covariance:
            for k in range(3, self.n_rbf + 1):
                rows, cols = np.indices((W.size, W.size))            
                row_vals = np.diag(rows, k=k)
                col_vals = np.diag(cols, k=k)
                corrcoef[row_vals, col_vals] = 1
                corrcoef[col_vals, row_vals] = 1
        
        corrcoef[self.n_rbf:, self.n_rbf:] = 0
            
        reg = corrcoef * np.atleast_2d(W).T * np.atleast_2d(W)
                
        def loss(θ):
            χ2 = np.sum((A @ θ - flux)**2 * ivar)
            #L2 = np.sum(Lambda * W * θ**2)
            L2 = np.sum(reg @ θ**2)
            return χ2 + L2
                
        def jacobian_loss(θ):
            d_chi2_d_theta = 2 * A.T @ ((A @ θ - flux) * ivar)
            d_L2_d_theta = 2 * reg @ θ #2 * Lambda * W * θ
            return d_chi2_d_theta + d_L2_d_theta 
        
        x_opt, nfeval, rc = op.fmin_tnc(
            loss,
            x0=x0,
            fprime=jacobian_loss,
            bounds=bounds.T,
            maxfun=10_000,
            disp=0,
            messages=0,
        )
        
        result_message = {
            -1: "Infeasible (lower bound > upper bound)",
            0: "Local minimum reached (|pg| ~= 0)",
            1: "Converged (|f_n-f_(n-1)| ~= 0)",
            2: "Converged (|x_n-x_(n-1)| ~= 0)",
            3: "Max. number of function evaluations reached",
            4: "Linear search failed",
            5: "All lower bounds are equal to the upper bounds",
            6: "Unable to progress",
            7: "User requested end of minimization",
        }.get(rc)
        
        if ax is None:
            fig, ax = plt.subplots()
            
        
        ax.plot(λ, flux, c='k')
        #ax.plot(λ, A @ p_init.x, c="tab:red")
        ax.plot(λ, A @ x_opt, c="tab:blue")
        
        central_mask = np.zeros(x0.size, dtype=bool)
        central_mask[self.n_rbf // 2 - central:self.n_rbf // 2 + central] = True
        central_mask[self.n_rbf:] = True
        ax.plot(λ, A[:, central_mask] @ x_opt[central_mask], c="tab:green")
        
        ax.set_ylim(0.5, 1.1)
        
        print(f"{loss(x0):.3e} -> {loss(x_opt):.3e}")
        
        # find contiguous places
        is_non_zero = (x_opt[:self.n_rbf] > 0)
        edges = 1 + np.where(np.diff(is_non_zero))[0]
        indices = np.sort(np.hstack([
            0, np.repeat(edges, 2), self.n_rbf 
        ]).reshape((-1, 2)))
        
        non_zero_indices = indices[0 if is_non_zero[0] else 1::2]
        for si, ei in non_zero_indices:
            assert is_non_zero[si:ei].all()
            assert si != ei
            v = np.sum(locations[si:ei] * x_opt[si:ei] / np.sum(x_opt[si:ei]))
            ew = 1000 * np.sum(x_opt[si:ei])
            if ew > 5:
                print(f"{si} {ei} {v:.3f} {ew:.3f}")
                ax.axvline(v, c="#666666", ls=":", zorder=-1)
                ax.axvspan(locations[si], locations[ei-1], color="#cccccc", alpha=0.5, zorder=-1)

        print()
        
        
    
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

    from tqdm import tqdm
    for j, λ_line in enumerate(tqdm(wls)):

        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True)

        #model = RBFSpectralLineModel(λ_line, continuum_basis=PolynomialBasis(0))
        #model.fit(spectrum, Lambda=1e6, regularize_covariance=True, ax=axes[0])
        
        model = RBFSpectralLineModel(λ_line, continuum_basis=PolynomialBasis(0))
        model.fit(spectrum, ax=axes[1], Λ=1e4)
    