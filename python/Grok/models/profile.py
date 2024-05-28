import numpy as np
from scipy import stats, special


class Profile:
    ...
        
class GaussianProfile(Profile):
    
    n_parameters = 3
    
    def __init__(self, λ, λ_tolerance=0.1):
        self.λ = λ
        self.λ_tolerance = λ_tolerance
        return None
    
    def __call__(self, λ, μ, σ, amplitude, *args, **kwargs):
        return 1 - amplitude * stats.norm.pdf(λ, loc=μ, scale=σ)
    
    
    def get_initial_guess(self, λ, flux, ivar):        
        return [self.λ, 0.05, 0.1]
                
    def get_bounds(self, σ_max=0.1, amplitude_max=1, **kwargs): #TODO: Put these to __init__ instead
        return [
            (self.λ - self.λ_tolerance, self.λ + self.λ_tolerance),
            (0, σ_max),
            (0, amplitude_max)
        ]
        
    def get_equivalent_width(self, θ):
        """Return the equivalent width of the fitted profile in milliAngstroms (mÅ)."""
        return 1000 * θ[2]

    def get_equivalent_width_uncertainty(self, θ, Σ):
        """Return the uncertainty on the equivalent width of the fitted profile in milliAngstroms (mÅ)."""
        return 1000 * Σ[2, 2]**0.5



class VoigtProfile(Profile):
    
    n_parameters = 4
    
    def __init__(self, λ, λ_tolerance=0.1):
        self.λ = λ
        self.λ_tolerance = λ_tolerance
        return None

    def get_initial_guess(self, λ, flux, ivar):
        return [self.λ, 0.1, 0, 0.2] # mu, sigma, gamma, amplitude
    
    def get_bounds(self, σ_max=0.1, amplitude_max=1, **kwargs):
        return [
            (self.λ - self.λ_tolerance, self.λ + self.λ_tolerance),
            (0, σ_max),    # sigma
            (0, np.inf),    # gamma
            (0, amplitude_max)          # amplitude
        ]
    
    def __call__(self, λ, μ, σ, γ, amplitude, *args, **kwargs):
        return 1 - amplitude * special.voigt_profile(λ - μ, σ, γ)
    
    def get_equivalent_width(self, θ, n_sigma=10, n_steps=1000, **kwargs):
        """Return the equivalent width of the fitted profile in milliAngstroms (mÅ)."""        
        μ, σ, *_ = θ        
        λ = np.linspace(μ - n_sigma*σ, μ + n_sigma*σ, n_steps)
        return 1000 * np.trapz(1 - self.__call__(λ, *θ), λ)
    
    def get_equivalent_width_uncertainty(self, θ, Σ, draws=1000, **kwargs):
        """Return the uncertainty on the equivalent width of the fitted profile in milliAngstroms (mÅ)."""        
        X = np.random.multivariate_normal(
            θ[:self.n_parameters], 
            Σ[:self.n_parameters, :self.n_parameters], 
            size=draws
        )
        p = np.array([self.get_equivalent_width(x, **kwargs) for x in X])
        return 0.5 * np.diff(np.percentile(p, [16, 84]))[0]
        
        

