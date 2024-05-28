import numpy as np


import numpy as np
import gzip
import warnings
import pickle
import os
from functools import cache
from typing import Optional, Tuple
from scipy import optimize as op
from sklearn.decomposition._nmf import _fit_coordinate_descent
from sklearn.exceptions import ConvergenceWarning
from scipy.signal.windows import tukey
from scipy.optimize import leastsq
from scipy import interpolate


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
    
    
    
    
SMALL = 1e-12

def expand_path(path):
    return os.path.abspath(os.path.expanduser(path))

@cache
def load_basis_vectors(path, P, pad=0):    
    full_path = expand_path(path)
    if full_path.lower().endswith(".gz"):
        with gzip.open(full_path, "rb") as fp:
            masked_basis_vectors = pickle.load(fp)
    else:        
        with open(full_path, "rb") as fp:
            masked_basis_vectors = pickle.load(fp)
    if pad > 0:
        basis_vectors = np.zeros((masked_basis_vectors.shape[0], P + 2 * pad))
        basis_vectors[:, pad:-pad] = masked_basis_vectors
    else:
        basis_vectors = masked_basis_vectors
    return basis_vectors

def fit_sines_and_cosines(wavelength, flux, ivar, deg: Optional[int] = 3, L: Optional[float] = None, A = None):
    L = L or np.ptp(wavelength)
    if A is None:
        A = design_matrix(wavelength, deg)
    MTM = A @ (ivar[:, None] * A.T)
    MTy = A @ (ivar * flux)
    try:
        theta = np.linalg.solve(MTM, MTy)
    except np.linalg.LinAlgError:
        if np.any(ivar > 0):
            raise
    continuum = A @ theta
    return (continuum, theta, A)

    
class BaseContinuumModel(object):

    def __init__(
        self,
        wavelength: np.array,
        basis_vectors: np.ndarray,
        deg: Optional[int] = 3,
    ):       
        self.wavelength = wavelength
        _check_wavelength_basis_vectors_shape(wavelength, basis_vectors)
        self.basis_vectors = basis_vectors
        self.deg = deg
        return None

    def fit_sines_and_cosines(self, wavelength, flux, ivar, A=None):
        if A is None:
            A = design_matrix(wavelength, self.deg)
        MTM = A @ (ivar[:, None] * A.T)
        MTy = A @ (ivar * flux)
        try:
            theta = np.linalg.solve(MTM, MTy)
        except np.linalg.LinAlgError:
            if np.any(ivar > 0):
                raise
        continuum = theta @ A
        return (continuum, theta, A)

    def _W_step(self, mean_rectified_flux, W, **kwargs):
        absorption = 1 - mean_rectified_flux
        use = np.ones(mean_rectified_flux.size, dtype=bool)
        use *= (
            np.isfinite(absorption) 
        &   (absorption >= 0) 
        &   (mean_rectified_flux > 0)
        )
        W_next, _, n_iter = _fit_coordinate_descent(
            absorption[use].reshape((1, -1)),
            W,
            self.basis_vectors[:, use],
            update_H=False,
            verbose=False,
            shuffle=True
        )        
        rectified_model_flux = 1 - (W_next @ self.basis_vectors)[0]
        return (W_next, rectified_model_flux, np.sum(use), n_iter)


    def _predict(self, theta, A_slice, C, P):
        return (1 - theta[:C] @ self.basis_vectors) * (A_slice @ theta[C:]).reshape((-1, P))

    
    def full_design_matrix(self, N, R=1):        
        C, P = self.basis_vectors.shape
        
        K = R * (2 * self.deg + 1)
        A = np.zeros((N * P, C + N * K), dtype=float)
        for i in range(N):
            A[i*P:(i+1)*P, :C] = self.basis_vectors.T
            A[i*P:(i+1)*P, C + i*K:C + (i+1)*K] = design_matrix(self.wavelength, self.deg).T
        return A

    def get_mask(self, ivar):
        N, P = np.atleast_2d(ivar).shape        
        use = np.zeros((N, P), dtype=bool)
        #use[:, np.hstack(self.region_masks)] = True
        use *= (ivar > 0)
        return ~use

    

    def _get_initial_guess_with_small_W(self, si, ei, flux, ivar, initial_theta=None, small=1e-12):
        with warnings.catch_warnings():        
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            N = 1
            C, P = self.basis_vectors.shape
            R = 1
            P = ei - si
            
            K = R * (2 * self.deg + 1)
            A = np.zeros((N * P, C + N * K), dtype=float)
            for i in range(N):
                A[i*P:(i+1)*P, :C] = self.basis_vectors[:, si:ei].T
                A[i*P:(i+1)*P, C + i*K:C + (i+1)*K] = design_matrix(self.wavelength[si:ei], self.deg).T


            Y = flux.flatten()
            use = ~self.get_mask(ivar).flatten()
            result = op.lsq_linear(
                A[use],
                Y[use],
                bounds=self.get_bounds(A.shape[1], [-np.inf, 0]),
            )            
            C, P = self.basis_vectors.shape       
            theta = result.x[C:]     
            if initial_theta is not None and np.allclose(theta, np.zeros_like(theta)):
                theta = initial_theta
            return np.hstack([small * np.ones(C), theta])
                
    
    def _get_initial_guess(self, si, ei, flux, continuum, initial_theta):
                
        absorption = np.clip(1 - flux / continuum, 0, 1)
        use = np.ones(absorption.size, dtype=bool)
        use *= (
            np.isfinite(absorption) 
        &   (absorption >= 0) 
        )
        W = SMALL * np.ones((1, self.basis_vectors.shape[0]))
        
        W_next, _, n_iter = _fit_coordinate_descent(
            absorption[use].reshape((1, -1)),
            W,
            self.basis_vectors[:, si:ei][:, use],
            update_H=False,
            verbose=False,
            shuffle=True
        )               
        return np.hstack([W_next.flatten(), initial_theta])


    



    def fit(self, wavelength, flux, ivar, v_rel=0.0):
        """
        Fit the continuum in every order simultaneously.
        """

        si, ei = self.wavelength.searchsorted(overlap(self.wavelength, wavelength))
        interp_wavelength = self.wavelength[si:ei]
        interp_flux = np.interp(interp_wavelength, wavelength, flux, left=0, right=0)
        interp_ivar = np.interp(interp_wavelength, wavelength, ivar, left=0, right=0)
        
        sigma = interp_ivar**-0.5

        C, P = self.basis_vectors.shape

        A = design_matrix(interp_wavelength, self.deg)
        (initial_continuum, initial_theta, _) = self.fit_sines_and_cosines(
            interp_wavelength, interp_flux, interp_ivar, 
            A=A
        )
        
        x0 = self._get_initial_guess(si, ei, interp_flux, initial_continuum, initial_theta)
        x0 = self._get_initial_guess_with_small_W(si, ei, interp_flux, interp_ivar, initial_theta=initial_theta)
        x0 = np.zeros_like(x0)
        x0[:C] = 10**-2

        basis_vectors_slice = self.basis_vectors[:, si:ei]
        
        def f(_, *theta):
            print(np.log10(theta[:C]))
            return (1 - theta[:C] @ basis_vectors_slice) * (A.T @ theta[C:]).flatten()
        
        p_opt, cov = op.curve_fit(
            f,
            None,
            interp_flux,
            p0=x0,
            sigma=sigma,
            bounds=self.get_bounds(C + 2 * self.deg + 1)
        )

        '''
        fig, ax = plt.subplots()
        ax.plot(interp_wavelength, interp_flux, c='k')
        ax.plot(interp_wavelength, f(None, *x0), c="tab:red")
        ax.plot(interp_wavelength, f(None, *x1), c="tab:blue")
        ax.plot(interp_wavelength, f(None, *p_opt), c="tab:green")
        '''

        rectified_model_flux = 1 - p_opt[:C] @ basis_vectors_slice
        # compute the continuum using the original design matrix

        # TODO: this is not quite right... change it all to use design matrices in one reference frame!!
        #continuum = np.interp(wavelength * (1 - initial_v_rel/3e5), interp_wavelength, p_opt[C:] @ A.T, left=0, right=0)
        
        interp_continuum = A.T @ p_opt[C:]
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(interp_wavelength, interp_flux, c='k')
        ax.plot(interp_wavelength, interp_continuum, c="tab:red")
        #init_continuum = A.T @ x0[C:]
        
        ax.plot(interp_wavelength, initial_continuum, c="tab:blue")
        
        #ax.plot(interp_wavelength, rectified_model_flux, c="tab:red")
        
        init_continuum = A.T @ x0[C:]
        init_rectified = 1 - x0[:C] @ basis_vectors_slice
        
        fig, ax = plt.subplots()
        ax.plot(interp_wavelength, rectified_model_flux, c="tab:red")
        ax.plot(interp_wavelength, init_rectified, c='tab:blue')
        
        
        
        raise a
        


    def fit_relative_velocity(self, wavelength, flux, ivar, x0=None, **kwargs):
        
        # TODO: check that wavelength, flux, etc are 1D arrays only
    
        A_obs = design_matrix(wavelength, self.deg)
        (initial_continuum, initial_theta, _) = self.fit_sines_and_cosines(wavelength, flux, ivar, A=A_obs)

        si, ei = np.searchsorted(self.wavelength, wavelength[[0, -1]])
        
        C, P = self.basis_vectors.shape
        v_rels, chi2s = (np.nan * np.ones(C), np.nan * np.ones(C))
        for i, basis_vector in enumerate(self.basis_vectors[:, si:ei]):
            v_rel, e_v_rel, chi2, corr, lags = measure_relative_velocity(
                wavelength,
                flux / initial_continuum,
                self.wavelength[si:ei],
                1 - basis_vector/np.max(basis_vector)
            )
            chi2s[i] = chi2
            v_rels[i] = v_rel
        
        initial_v_rel = v_rels[np.argmin(chi2s)]
        
        # interpolate data to model (this is bad practice, but all of this is approximate)
        si, ei = self.wavelength.searchsorted(overlap(self.wavelength, wavelength))
        interp_wavelength = self.wavelength[si:ei]
        interp_flux = np.interp(interp_wavelength, wavelength * (1 - initial_v_rel/3e5), flux, left=0, right=0)
        interp_ivar = np.interp(interp_wavelength, wavelength * (1 - initial_v_rel/3e5), ivar, left=0, right=0)
        interp_continuum = np.interp(interp_wavelength, wavelength * (1 - initial_v_rel/3e5), initial_continuum, left=0, right=0)
        
        sigma = interp_ivar**-0.5
        
        x0 = self._get_initial_guess(si, ei, interp_flux, interp_continuum, initial_theta)
        
        A = design_matrix(interp_wavelength, self.deg).T

        basis_vectors_slice = self.basis_vectors[:, si:ei]
        
        def f(_, *theta):            
            return (1 - theta[:C] @ basis_vectors_slice) * (A @ theta[C:]).flatten()
        
        p_opt, cov = op.curve_fit(
            f,
            None,
            interp_flux,
            p0=x0,
            sigma=sigma,
            bounds=self.get_bounds(C + 2 * self.deg + 1)
        )

        '''
        fig, ax = plt.subplots()
        ax.plot(interp_wavelength, interp_flux, c='k')
        ax.plot(interp_wavelength, f(None, *x0), c="tab:red")
        ax.plot(interp_wavelength, f(None, *x1), c="tab:blue")
        ax.plot(interp_wavelength, f(None, *p_opt), c="tab:green")
        '''

        rectified_model_flux = 1 - p_opt[:C] @ basis_vectors_slice
        # compute the continuum using the original design matrix
        # TODO: this is not quite right... change it all to use design matrices in one reference frame!!
        continuum = np.interp(wavelength * (1 - initial_v_rel/3e5), interp_wavelength, p_opt[C:] @ A.T, left=0, right=0)
        
        
        # measure the radial velocity again
        v_rel, e_v_rel, chi2, corr, lags = measure_relative_velocity(
            wavelength,
            flux / continuum,
            interp_wavelength,
            rectified_model_flux
        )
        
        return (v_rel, continuum, interp_wavelength, rectified_model_flux, p_opt)


    def get_bounds(self, N, component_bounds=(0, +np.inf)):
        C, P = self.basis_vectors.shape          

        return np.vstack([
            np.tile(component_bounds, C).reshape((C, 2)),
            np.tile([-np.inf, +np.inf], N - C).reshape((-1, 2))
        ]).T            


def overlap(A, B):
    # get range of values in A and B that overlap
    return (np.max([A[0], B[0]]), np.min([A[-1], B[-1]]))






def _check_wavelength_basis_vectors_shape(wavelength, basis_vectors):
    P = wavelength.size
    assert wavelength.ndim == 1, "wavelength must be a one-dimensional array." 
    C, P2 = basis_vectors.shape
    assert P == P2, "`basis_vectors` should have shape (C, P) where P is the size of `wavelength`"


def design_matrix(wavelength: np.array, deg: int) -> np.array:
    L = 2 * np.ptp(wavelength)
    scale = 2 * (np.pi / L)
    return np.vstack(
        [
            np.ones_like(wavelength).reshape((1, -1)),
            np.array(
                [
                    [np.cos(o * scale * wavelength), np.sin(o * scale * wavelength)]
                    for o in range(1, deg + 1)
                ]
            ).reshape((2 * deg, wavelength.size)),
        ]
    )    
    
    
# self calib stuff
import numpy as np
from typing import Union, Tuple, Optional
from functools import cached_property

from scipy import optimize as op
from sklearn.decomposition._nmf import _initialize_nmf

EPSILON = np.finfo(np.float32).eps

class BaseSelfCalibratedContinuum(object):

    def __init__(
        self,
        wavelength: np.array,
        deg: int,
        L: Union[float, int],
        regions: Optional[Tuple[Tuple[float, float]]] = None,
        tol: float = 1e-10,
        n_components: int = 1,
        outer_iter: int = 1000,
        inner_W_iter: int = 1,
        inner_WH_iter: int = 1,
        init: Optional[str] = None,  
        initial_continuum_scalar = 1,   
        alpha_W=0,
        alpha_H=0,
        l1_ratio=0,
        **kwargs
    ):    
        self.L = L or wavelength.size
        self.deg = deg
        self.tol = tol
        trim = 20
        self.alpha_W = alpha_W
        self.alpha_H = alpha_H
        self.l1_ratio = l1_ratio
        self.regions = regions or (wavelength[[trim, -trim]], )
        self.wavelength = wavelength
        self.n_components = n_components
        self.outer_iter = outer_iter
        self.inner_W_iter = inner_W_iter
        self.inner_WH_iter = inner_WH_iter
        self.init = init
        self.initial_continuum_scalar = initial_continuum_scalar
        return None
    
    @cached_property
    def region_masks(self):
        return _region_masks(self.wavelength, self.regions)
    
    @cached_property
    def continuum_design_matrix(self):
        return _continuum_design_matrix(self.wavelength, self.deg, self.L, self.regions, self.region_masks)
    
    @cached_property
    def n_parameters_per_region(self):
        #return self.deg + 1
        return 2 * self.deg + 1
    
    def _fit_continuum(self, flux, ivar):
        N, P = flux.shape
        theta = np.zeros((N, len(self.regions), self.n_parameters_per_region))
        continuum = np.nan * np.ones_like(flux)
        for i in range(N):            
            for j, mask in enumerate(self.region_masks):
                sj, ej = (j * self.n_parameters_per_region, (j + 1) * self.n_parameters_per_region)
                A = self.continuum_design_matrix[mask, sj:ej]
                MTM = A.T @ (ivar[i, mask][:, None] * A)
                MTy = A.T @ (ivar[i, mask] * flux[i, mask])
                try:
                    theta[i, j] = np.linalg.solve(MTM, MTy)
                except np.linalg.LinAlgError:
                    if np.any(ivar[i, mask] > 0):
                        print(f"Continuum warning on {i}, {j}")
                        continue
                continuum[i, mask] = A @ theta[i, j]        

        return (theta, continuum)
    
    '''
    def _fit_continuum(self, flux, ivar):
        N, P = flux.shape
        theta = np.zeros((N, len(self.regions), self.n_parameters_per_region))
        continuum = np.nan * np.ones_like(flux)
        for i in range(N):
            for j, mask in enumerate(self.region_masks):
                A = np.vander(self.wavelength[mask], self.deg + 1)
                C_inv = np.diag(ivar[i, mask])
                MTM = A.T @ C_inv @ A
                MTy = A.T @ C_inv @ flux[i, mask]
                theta[i, j] = np.linalg.solve(MTM, MTy)
                continuum[i, mask] = A @ theta[i, j]
                
        return (theta, continuum)
    '''
    
    def _fit_W_and_theta_by_optimization(self, flux, ivar, H, W_init, theta_init):
        N, P = flux.shape
        theta = np.copy(theta_init)
        W = np.copy(W_init)
        continuum = np.zeros_like(flux)
        
        use = ~self.get_mask(ivar)
        sigma = ivar**-0.5
        chi2, dof = (0, 0)
        
        for i in range(N):
            
            mask = use[i]
            
            def f(_, *p):
                return (
                    (1 - p[:self.n_components] @ H) 
                *   (self.continuum_design_matrix @ p[self.n_components:])
                )[mask]

            p0 = np.hstack([W_init[i], theta_init[i].ravel()])
            
            try:
                p_opt, cov = op.curve_fit(
                    f,
                    None,
                    flux[i][mask],
                    p0=p0,
                    sigma=sigma[i][mask],
                    bounds=self.get_bounds(1)
                )
            except KeyboardInterrupt:
                raise
            except:
                print(f"failed on {i}")
                p_opt = p0
            else:
                W[i] = p_opt[:self.n_components]
                theta[i] = p_opt[self.n_components:].reshape((len(self.regions), self.n_parameters_per_region))
            
            finally:
                continuum[i] = self.continuum_design_matrix @ p_opt[self.n_components:]
                chi2 += np.sum((f(None, *p_opt) - flux[i][mask])**2 * ivar[i][mask])
                dof += np.sum(mask)
                
        return (W, theta, continuum, chi2, dof)


    def _fit_W_and_theta_by_step(self, flux, ivar, X, V, W, H, inner_W_iter=None, l1_reg_H=0, l2_reg_H=0, l1_reg_W=0, l2_reg_W=0):
    
        N, P = X.shape        
        W = np.copy(W)
        theta = np.zeros((N, len(self.regions), self.n_parameters_per_region))
        continuum = np.zeros_like(flux)
                
        chi2, dof = (0, 0)        
        for i in range(N):
            
            # fit W
            for _ in range(inner_W_iter or self.inner_W_iter):
                _multiplicative_update(X, V, W, H, update_H=False, update_W=True, l1_reg_H=l1_reg_H, l2_reg_H=l2_reg_H, l1_reg_W=l1_reg_W, l2_reg_W=l2_reg_W)
                #W[i] *= ((X[[i]] * V[[i]] @ H.T) / ((V[[i]] * (W[[i]] @ H)) @ H.T))[0]
            
            # fit theta
            model_flux = 1 - (W[[i]] @ H)
            continuum_flux = flux[i] / model_flux
            continuum_ivar = ivar[i] * model_flux ** 2
            
            theta_, continuum_ = self._fit_continuum(continuum_flux, continuum_ivar)
            continuum[i] = continuum_
            theta[i] = theta_
            '''
            for j, mask in enumerate(self.region_masks):                
                sj, ej = (j * self.n_parameters_per_region, (j + 1) * self.n_parameters_per_region)
                A = self.continuum_design_matrix[mask, sj:ej]
                MTM = A.T @ (continuum_ivar[0, mask][:, None] * A)
                MTy = A.T @ (continuum_ivar[0, mask] * continuum_flux[0, mask])
                try:
                    theta[i, j] = np.linalg.solve(MTM, MTy)
                except np.linalg.LinAlgError:
                    if np.any(continuum_ivar[0, mask] > 0):
                        print(f"Continuum warning on {i}, {j}")
                        continue
                continuum[i, mask] = A @ theta[i, j]        
            '''
            
            diff = (flux[i] - continuum[i] * model_flux[0])**2 * ivar[i]
            finite = np.isfinite(diff)
            chi2 += np.sum(diff[finite])
            dof += np.sum(finite)
            
        #assert dof > 0
        
        return (W, theta, continuum, chi2, dof)


    def get_mask(self, ivar):
        N, P = np.atleast_2d(ivar).shape        
        use = np.zeros((N, P), dtype=bool)
        use[:, np.hstack(self.region_masks)] = True
        use *= (ivar > 0)
        return ~use            

    def get_bounds(self, N, component_bounds=(0, +np.inf)):
        C = self.n_components
        A = N * len(self.regions) * (self.n_parameters_per_region)
        return np.vstack([
            np.tile(component_bounds, C).reshape((C, 2)),
            np.tile([-np.inf, +np.inf], A).reshape((-1, 2))
        ]).T            


    def fit(self, flux, ivar=None, initial_H=None):
        
        flux, ivar = _check_and_reshape_flux_ivar(flux, ivar, self.wavelength)
        
        l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H = _compute_regularization(flux, self.alpha_W, self.alpha_H, self.l1_ratio)

        # (1) fit theta (no eigenvectors)
        # (2) fit eigenvectors given theta
        # (3) fit W given H, theta
        # (4) fit theta given W, H
        # (5) go to step 2
        
        # allow for an initial H
        
        theta, continuum = self._fit_continuum(flux, ivar)
        
        continuum *= self.initial_continuum_scalar
        
        X, V = _get_XV(flux, ivar, continuum)
        if initial_H is None:
            W, H = _initialize_nmf(X, self.n_components, init=self.init, eps=1e-6, random_state=None)
            W = np.atleast_2d([1]).astype(W.dtype)
            H = np.clip(1 - np.atleast_2d(flux / continuum), EPSILON, None).astype(H.dtype)
            H[~np.isfinite(H)] = EPSILON
        else:
            print("using initial H")
            H = np.copy(initial_H)
            W = 1e-12 * np.ones((flux.shape[0], self.n_components))
            
            W, theta, continuum, chi2, dof = self._fit_W_and_theta_by_step(
                flux, ivar, X, V, W, H, inner_W_iter=1,
                l1_reg_H=l1_reg_H, l2_reg_H=l2_reg_H, l1_reg_W=l1_reg_W, l2_reg_W=l2_reg_W
                )
            X, V = _get_XV(flux, ivar, continuum)
            print(f"initial chi2: {chi2/dof}")
        
        last_rchi2 = None
        try:
            #with (total=self.outer_iter) as pb:
            for outer in range(self.outer_iter):    
                for _ in range(self.inner_WH_iter): # sub-iters
                    _multiplicative_update(X, V, W, H, update_W=True, update_H=True, l1_reg_H=l1_reg_H, l2_reg_H=l2_reg_H, l1_reg_W=l1_reg_W, l2_reg_W=l2_reg_W)
                
                # For each spectrum, fit W + theta given H
                W, theta, continuum, chi2, dof = self._fit_W_and_theta_by_step(flux, ivar, X, V, W, H, l1_reg_H=l1_reg_H, l2_reg_H=l2_reg_H, l1_reg_W=l1_reg_W, l2_reg_W=l2_reg_W)
                
                X, V = _get_XV(flux, ivar, continuum)
                
                #pb.set_description(f"chi2/dof={chi2/dof:.2e}")
                #pb.update()
                rchi2 = chi2 / dof
                if last_rchi2 is None:
                    last_rchi2 = rchi2
                else:
                    diff = rchi2 - last_rchi2
                    last_rchi2 = rchi2
                    print(outer, diff)
                    
                    if (diff <= 0 and np.abs(diff) < self.tol) or (diff > 0):
                        break
                    
        except KeyboardInterrupt:
            print("Detected KeyboardInterrupt: finishing the fit early")
            None

        return (W, H, theta, continuum, chi2, dof)
            

def _compute_regularization(X, alpha_W, alpha_H, l1_ratio):
    n_samples, n_features = X.shape
    l1_reg_W = n_features * alpha_W * l1_ratio
    l1_reg_H = n_samples * alpha_H * l1_ratio
    l2_reg_W = n_features * alpha_W * (1.0 - l1_ratio)
    l2_reg_H = n_samples * alpha_H * (1.0 - l1_ratio)   
    return (l1_reg_W, l1_reg_H, l2_reg_W, l2_reg_H)
        
def _multiplicative_update(X, V, W, H, update_H=True, update_W=True, l1_reg_H=0, l2_reg_H=0, l1_reg_W=0, l2_reg_W=0):
    WH = W @ H
    if update_H:
        numerator = ((V.T * X.T) @ W).T
        denominator = ((V.T * WH.T) @ W).T
        if l1_reg_H > 0:
            denominator += l1_reg_H
        if l2_reg_H > 0:
            denominator += l2_reg_H * H
        denominator[denominator == 0] = EPSILON
        H *= numerator / denominator
    
    if update_W:
        numerator = (X * V) @ H.T
        denominator = (V * WH) @ H.T
        if l1_reg_W > 0:
            denominator += l1_reg_W
        if l2_reg_W > 0:
            denominator += l2_reg_W * W
        W *= numerator/denominator
    return None


def _get_XV(flux, ivar, continuum):
    X = 1 - flux / continuum
    V = ivar * continuum**2
    
    is_bad_pixel = (
        (V == 0)
    |   (~np.isfinite(V))
    |   (X < 0)
    |   (~np.isfinite(X))
    )
    X[is_bad_pixel] = V[is_bad_pixel] = 0
    return (X, V)


def _design_matrix(wavelength: np.array, deg: int, L: float) -> np.array:
    scale = 2 * (np.pi / L)
    return np.vstack(
        [
            np.ones_like(wavelength).reshape((1, -1)),
            np.array(
                [
                    [np.cos(o * scale * wavelength), np.sin(o * scale * wavelength)]
                    for o in range(1, deg + 1)
                ]
            ).reshape((2 * deg, wavelength.size)),
        ]
    )


def _continuum_design_matrix(wavelength, deg, L, regions, region_masks):
    n_parameters_per_region = 2 * deg + 1    
    A = np.zeros((wavelength.size, len(regions) * n_parameters_per_region), dtype=float)
    for i, mask in enumerate(region_masks):
        si = i * n_parameters_per_region
        ei = (i + 1) * n_parameters_per_region
        A[mask, si:ei] = _design_matrix(wavelength[mask], deg, L).T        
    return A


def _region_masks(wavelength, regions):
    slices = []
    for region in regions:
        si, ei = wavelength.searchsorted(region)
        slices.append(np.arange(si, ei + 1, dtype=int))
    return slices


def _check_and_reshape_flux_ivar(flux, ivar, wavelength=None):
    if ivar is None:
        ivar = np.ones_like(flux)
    flux, ivar = (np.atleast_2d(flux), np.atleast_2d(ivar))
    N1, P1 = flux.shape
    N2, P2 = ivar.shape

    assert (N1 == N2) and (P1 == P2), "`flux` and `ivar` do not have the same shape"
    if wavelength is not None:
        P = len(wavelength)
        assert (P == P1), f"Number of pixels in flux does not match wavelength array ({P} != {P1})"

    bad_pixel = (
        (~np.isfinite(flux))
    |   (~np.isfinite(ivar))
    |   (flux <= 0)
    )
    flux[bad_pixel] = 0
    ivar[bad_pixel] = 0
    return (flux, ivar)


def fit_nmf_sinusoids(wavelength, flux, ivar, deg, L, **kwargs):
    
    
    
    model = BaseSelfCalibratedContinuum(wavelength=wavelength, deg=deg, L=L, **kwargs)
    
    # ignore 0 flux entries
    '''
    flux = flux.copy()
    ivar = ivar.copy()

    bad_pixel = (flux <= 0)    
    bad_pixel[(3940 >= wavelength) * (wavelength >= 3930)] = True
    bad_pixel[(3980 >= wavelength) * (wavelength >= 3960)] = True
    bad_pixel[(7700 >= wavelength) * (wavelength >= 7575)] = True
    
    ivar[bad_pixel] = 0
    ''' 
    return model.fit(flux, ivar)
    
    
    

if __name__ == "__main__":
    
    index = 30
    model_wavelength = 10 * (10**(2.57671464 + np.arange(167283) * 2.25855074e-06))
    basis_vectors = load_basis_vectors("/Users/andycasey/research/Grok/H.pkl.gz", 0)

    model = BaseContinuumModel(model_wavelength, basis_vectors, deg=2)

    from spectrum import Spectrum1D
    
    spectra = Spectrum1D.read("/Users/andycasey/Downloads/hd122563blue_multi.fits")
    #spectra = Spectrum1D.read("/Users/andycasey/Downloads/j174239-133332blue_multi.fits")

    v_rel = -42
    
    def air_to_vacuum(lambdas):
        """
        Convert air wavelengths to vacuum.

        As per: https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
        
        :param lambdas:
            The air wavelengths.
        """
        
        values = lambdas#.to("Angstrom").value
        s = 10**4 / values
        n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
        return values * n #* unit
        
    import gzip
    import pickle
    
    with gzip.open("../../H.pkl.gz", "rb") as fp:
        foo = pickle.load(fp)
        
            
    wl, flux, ivar = ([], [], [])
    for index in (21, ):
        wl.extend(air_to_vacuum(spectra[index].wavelength * (1 - v_rel/299792.458)))
        flux.extend(spectra[index].flux)
        ivar.extend(spectra[index].ivar)
    #wl, flux, ivar = spectra[index].wavelength, spectra[index].flux, spectra[index].ivar
    
    #wl = wl * (1 - v_rel/299792.458)
    '''
    idx = np.argsort(wl)
    wl, flux, ivar = tuple(map(np.array, (wl, flux, ivar)))
    wl = wl[idx]
    flux = flux[idx]
    ivar = ivar[idx]
    '''
    
    model = BaseContinuumModel(model_wavelength, foo, 5)
    
    foo = model.fit(wl, flux, ivar)
    
    raise a
    
    
    #model.fit(wl, flux, ivar)
    
    #self_model = BaseSelfCalibratedContinuum(wavelength=wl, deg=3, L=wl.size)
    #(W, H, theta, continuum, chi2, dof) = self_model.fit(flux, ivar)
        

    
    #interp_basis_vectors = np.array([np.interp(model_wavelength )])
    
    
    
    import matplotlib.pyplot as plt
    H_reg = 1e3
    W_reg = 1e3
    (W, H, theta, continuum, chi2, dof) = fit_nmf_sinusoids(
        wl, 
        flux, 
        ivar, 
        5, 
        L=300,
        #initial_H=
        l1_reg_H=H_reg,
        l2_reg_H=H_reg,
        l1_reg_W=W_reg,
        l2_reg_W=W_reg,
        inner_WH_iter=10,
        n_components=1
)
    
    fig, ax = plt.subplots()
    #ax.plot(wl, flux / continuum[0], c='k')
    #ax.plot(wl, 1 - (W @ H)[0], c="tab:red")
    ax.plot(wl, flux, c='k')
    ax.plot(wl, (1 - (W @ H)[0]) * continuum[0], c='tab:red')
