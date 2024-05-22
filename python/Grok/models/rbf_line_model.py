
import numpy as np
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
        n_rbf_basis: int = 101,
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
        if (n_rbf_basis % 2) == 0:
            raise ValueError("n_rbf_basis must be odd")
        self.n_rbf_basis = n_rbf_basis
        return None
    
    def _prepare_fit(self, spectrum):
        # Get the data.        
        si, ei = spectrum.λ.searchsorted([self.λ - self.window, self.λ + self.window]) + [0, 1]
        λ, flux, ivar = (spectrum.λ[si:ei], spectrum.flux[si:ei], spectrum.ivar[si:ei])
    
        # Build design matrices.    
        locations = np.hstack([
            np.linspace(λ[0], self.λ, self.n_rbf_basis // 2 + 1)[:-1],
            np.linspace(self.λ, λ[-1], self.n_rbf_basis // 2 + 1)            
        ])
        
        rb = RadialBasis(locations, self.sigma)
        A_rb = rb.get_design_matrix(λ)
        bounds = [(0, np.inf)] * self.n_rbf_basis
        if self.continuum_basis is not None:
            A_c = self.continuum_basis.get_design_matrix(λ)
            bounds.extend(self.continuum_basis.get_bounds())
        else:
            A_c = np.empty((λ.size, 0))
        
        A = np.hstack([A_rb, A_c])
        return (λ, flux, ivar, locations, A_rb, A_c, A, np.atleast_2d(bounds).T)
                
                
    def _fit_least_squares(self, A, flux, ivar, Lambda, full_output=False):
        # Fit the model.
        ATCinv = A.T * ivar
        ATCinvA = ATCinv @ A
        ATCinvY = ATCinv @ flux
        R = Lambda * np.ones(A.shape[1])
        R[self.n_rbf_basis:] = 0 # do not regularize the continuum coefficients
        central = 3
        R[self.n_rbf_basis // 2 - central:self.n_rbf_basis // 2 + central] = 0 # do not regularize the central few bases 
        R = np.diag(R)

        θ = np.linalg.solve(ATCinvA + R, ATCinvY)
        
        if full_output:
            return (θ, ATCinv, ATCinvA, ATCinvY)
        return θ
            
        
    
    def fit(
        self,
        spectrum,
        op_kwds: Optional[dict] = None,
        Lambda: Optional[float] = 1e3,
        **kwargs
    ):
        λ, flux, ivar, locations, A_rb, A_c, A, bounds = self._prepare_fit(spectrum)

        kwds = dict(max_iter=10_000)
        kwds.update(op_kwds or {})
        
        # Let's get a bounded unregularized estimate that does not account for uncertainties.
        p_init = op.lsq_linear(A, flux, bounds, **kwds)        

        R = Lambda * np.ones_like(p_init.x)
        R[self.n_rbf_basis:] = 0
        central = 3
        R[self.n_rbf_basis // 2 - central:self.n_rbf_basis // 2 + central] = 0        
        
        def loss(θ):
            χ2 = np.sum((A @ θ - flux)**2 * ivar)
            #L1 = np.sum(np.abs(R * θ))
            L2 = np.sum(R * θ**2)
            return χ2 + L2
                
        def jacobian_loss(θ):
            """
            Calculates the Jacobian of the combined loss function with respect to θ.

            Args:
                θ: The parameter vector (NumPy array).
                A: The design matrix (NumPy array).
                flux: The observed flux values (NumPy array).
                ivar: The inverse variance weights (NumPy array).
                Lambda: The L1 regularization strength (scalar).

            Returns:
                The Jacobian matrix (NumPy array) with shape (len(θ),).
            """

            # Chi-squared term (χ²)
            d_chi2_d_theta = 2 * A.T @ ((A @ θ - flux) * ivar)

            # L1 term
            #d_L1_d_theta = R * np.sign(θ)
            d_L2_d_theta = 2 * R * θ

            # Combined Jacobian
            jacobian = d_chi2_d_theta + d_L2_d_theta #+ d_L1_d_theta

            return jacobian        
        
                
        ftol = gtol = np.finfo(float).eps
        result = self._fit_least_squares(A, flux, ivar, Lambda)
        
        foo = op.check_grad(loss, jacobian_loss, result, epsilon=1e-3)
        print(f"foo:{foo:.3e}")

        
        xopt, fopt, meta = op.fmin_l_bfgs_b(
            loss,
            fprime=jacobian_loss,
            x0=p_init.x, 
            bounds=bounds.T,
            #factr=1e3,
            iprint=True
        )
        print(meta)
        
        '''
        Y = flux
        Cinv = np.diag(ivar)
        W = R
        L = bounds[0]
        U = bounds[0]
        L[~np.isfinite(L)] = -1e16
        U[~np.isfinite(U)] = 1e16
        
        def projected_gradient_descent(A, X_init, learning_rate=1e-8, max_iterations=1000, tolerance=1e-6):
            """
            Solves the constrained linear system Y = AX with weighted regularization
            using projected gradient descent.

            Args:
                A: Design matrix (n_observations x n_variables)
                Y: Observation vector (n_observations)
                C_inv: Inverse covariance matrix (n_observations x n_observations)
                W: Diagonal matrix of regularization weights (n_variables x n_variables)
                L: Lower bounds on variables (n_variables)
                U: Upper bounds on variables (n_variables)
                learning_rate: Step size for gradient descent
                max_iterations: Maximum number of iterations
                tolerance: Convergence tolerance

            Returns:
                The optimal solution vector X (n_variables)
            """

            n_variables = A.shape[1]
            X = np.zeros(n_variables) + X_init  # Initial guess (within constraints)

            for _ in range(max_iterations):
                # Calculate gradient
                gradient = -A.T @ Cinv @ (Y - A @ X) + W @ X

                # Update X
                X_new = X - learning_rate * gradient

                # Project onto feasible region
                X_new = np.clip(X_new, L, U)
                
                assert np.all(np.isfinite(X_new))

                # Check for convergence
                if np.linalg.norm(X_new - X) < tolerance:
                    print(f"reached after {_} iters")
                    raise a
                    break

                X = X_new

            return X
        
        X_opt = projected_gradient_descent(A, p_init.x)
        '''
        
        X_opt, nfeval, rc = op.fmin_tnc(
            loss,
            x0=p_init.x,
            fprime=jacobian_loss,
            bounds=bounds.T,
        )

        '''
        result = op.least_squares(
            f,
            x0=p_init.x,
            bounds=bounds,
            verbose=2
        )
        '''
        
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(λ, flux, c='k')
        ax.plot(λ, A @ p_init.x, c="tab:red")
        #ax.plot(λ, A @ result, c="tab:blue")
        #ax.plot(λ, A @ xopt, c="tab:green")
        ax.plot(λ, A @ X_opt, c="tab:orange")
        
        '''
        for loc, v in zip(locations, p_init.x[:self.n_rbf_basis]):
            if v > 5e-3: # 5 mA
                ax.axvline(loc, c="#666666", ls=":")
        '''
        ax.set_ylim(0, 1.2)
        
        print(loss(p_init.x))
        #print(loss(xopt))
        print(loss(X_opt))
        raise a
    
        
        
    
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
    
    model = RBFSpectralLineModel(5068.7655, continuum_basis=PolynomialBasis(1))
    model.fit(spectrum, Lambda=1e10)
