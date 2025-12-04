"""
Fiber Modes: Core numerical solver for step-index optical fiber modes.

This module solves the LP₀₁ fundamental mode dispersion relation and computes
the mode profiles. It forms the foundation for perturbation analysis.

Theory Reference:
- Snyder & Love (1983): "Optical Waveguide Theory"
- Section 3-4 of PROJECT.md: Eigenmode theory and Sturm-Liouville formulation
"""

import numpy as np
from scipy.special import jv, kv, jvp, kvp  # Bessel functions and derivatives
from scipy.optimize import brentq, root_scalar
import warnings


class StepIndexFiber:
    """
    Represents a step-index optical fiber and solves for its guided modes.
    
    Parameters
    ----------
    n1 : float
        Core refractive index
    n2 : float
        Cladding refractive index
    a : float
        Core radius (in micrometers)
    wavelength : float
        Free-space wavelength λ₀ (in micrometers)
    """
    
    def __init__(self, n1, n2, a, wavelength):
        """Initialize fiber parameters."""
        if n1 <= n2:
            raise ValueError("Core index n1 must be > cladding index n2")
        if a <= 0 or wavelength <= 0:
            raise ValueError("Core radius and wavelength must be positive")
        
        self.n1 = n1
        self.n2 = n2
        self.a = a  # Core radius in µm
        self.wavelength = wavelength  # in µm
        
        # Precomputed convenience
        self._k0 = 2 * np.pi / wavelength  # Wave number (µm⁻¹)
        
        # Normalized frequency (V-number)
        self.V = self._k0 * a * np.sqrt(n1**2 - n2**2)
        
        # Will be computed on demand
        self._u = None
        self._w = None
        self._beta = None
        self._n_eff = None
    
    @property
    def k0(self):
        """Free-space wave number 2π/λ₀ (µm⁻¹)."""
        return self._k0
    
    def _dispersion_relation_lp01(self, u):
        """
        LP₀₁ fundamental mode dispersion relation.
        
        The eigenvalue equation is:
            J₀(u) / (u J₁(u)) + K₀(w) / (w K₁(w)) = 0
        
        We reformulate as:
            J₀(u) · w · K₁(w) + K₀(w) · u · J₁(u) = 0
        
        where w = √(V² - u²)
        
        Parameters
        ----------
        u : float
            Normalized radial eigenvalue (core region parameter)
        
        Returns
        -------
        float
            Residual of dispersion relation (should be zero at solution)
        """
        if u <= 0 or u >= self.V:
            return np.inf  # Unphysical region
        
        w_sq = self.V**2 - u**2
        if w_sq <= 0:
            return np.inf
        
        w = np.sqrt(w_sq)
        
        # Avoid singularities
        if w == 0 or u == 0:
            return np.inf
        
        try:
            # Bessel functions
            J0_u = jv(0, u)
            J1_u = jv(1, u)
            K0_w = kv(0, w)
            K1_w = kv(1, w)
            
            # Check for valid values
            if not (np.isfinite(J0_u) and np.isfinite(J1_u) and 
                    np.isfinite(K0_w) and np.isfinite(K1_w)):
                return np.inf
            
            # Avoid divide by zero
            if abs(J1_u) < 1e-15 or abs(K1_w) < 1e-15:
                return np.inf
            
            # Dispersion equation: J₀(u)/(u J₁(u)) + K₀(w)/(w K₁(w)) = 0
            # Reformulated to avoid division by small numbers:
            residual = J0_u * w * K1_w + K0_w * u * J1_u
            return residual
        except Exception as e:
            warnings.warn(f"Bessel function evaluation failed at u={u}: {e}")
            return np.inf
    
    def solve_dispersion(self, initial_guess=None, verbose=False):
        """
        Solve the LP₀₁ dispersion relation to find the eigenvalue u.
        
        Parameters
        ----------
        initial_guess : float, optional
            Initial guess for u. If None, uses V/2.
        verbose : bool, optional
            Print debug information.
        
        Returns
        -------
        dict
            Solution dictionary containing:
            - u : normalized core parameter
            - w : normalized cladding parameter
            - beta : propagation constant β (µm⁻¹)
            - n_eff : effective index n_eff = β/k₀
            - V : normalized frequency
        """
        if self.V < 0.01:
            raise ValueError("V-number too small for guided modes (V > 0)")
        
        # Initial guess
        if initial_guess is None:
            u_guess = 0.5 * self.V
        else:
            u_guess = initial_guess
        
        # Ensure u is in valid range [0, V)
        u_guess = np.clip(u_guess, 1e-3, self.V - 1e-3)
        
        # Find valid bracket by scanning
        u_brackets = np.linspace(0.01, self.V - 0.01, 200)
        f_brackets = []
        
        for u in u_brackets:
            f = self._dispersion_relation_lp01(u)
            f_brackets.append(f)
        
        f_brackets = np.array(f_brackets)
        
        # Find sign changes - be more permissive about what counts as finite
        u_lower, u_upper = None, None
        for i in range(len(f_brackets) - 1):
            f1 = f_brackets[i]
            f2 = f_brackets[i+1]
            
            # Both must be finite (not inf, not nan)
            if not (np.isinf(f1) or np.isnan(f1) or np.isinf(f2) or np.isnan(f2)):
                # Check for sign change
                if f1 * f2 < 0:
                    u_lower = u_brackets[i]
                    u_upper = u_brackets[i+1]
                    if verbose:
                        print(f"Found bracket: u ∈ [{u_lower:.6f}, {u_upper:.6f}], f ∈ [{f1:.6e}, {f2:.6e}]")
                    break
        
        if u_lower is None or u_upper is None:
            raise ValueError(f"No root found in bracket for V={self.V:.4f}")
        
        try:
            # Use Brent's method for robust root finding
            u_solution = brentq(
                self._dispersion_relation_lp01,
                u_lower,
                u_upper,
                xtol=1e-8,
                rtol=1e-10
            )
        except ValueError as e:
            # If Brent fails, try with bisect
            warnings.warn(f"Brent's method failed at [{u_lower}, {u_upper}]: {e}. Trying bisect.")
            try:
                result = root_scalar(
                    self._dispersion_relation_lp01,
                    bracket=[u_lower, u_upper],
                    method='bisect'
                )
                u_solution = result.root
            except Exception as e2:
                raise ValueError(f"Root finding failed: {e2}")
        
        # Compute w from V-parameter relation
        w_solution = np.sqrt(self.V**2 - u_solution**2)
        
        # Compute propagation constant β (in µm⁻¹)
        # β² = (k₀ n₁)² - u²/a²
        beta_solution = np.sqrt((self.k0 * self.n1)**2 - (u_solution / self.a)**2)
        
        # Effective index
        n_eff_solution = beta_solution / self.k0
        
        # Store results
        self._u = u_solution
        self._w = w_solution
        self._beta = beta_solution
        self._n_eff = n_eff_solution
        
        if verbose:
            residual = self._dispersion_relation_lp01(u_solution)
            print(f"✓ Dispersion equation solved:")
            print(f"  V-number: {self.V:.6f}")
            print(f"  u = {u_solution:.8f}")
            print(f"  w = {w_solution:.8f}")
            print(f"  Residual: {residual:.2e}")
            print(f"  β = {beta_solution:.8f} µm⁻¹")
            print(f"  n_eff = {n_eff_solution:.8f}")
        
        return {
            'u': u_solution,
            'w': w_solution,
            'beta': beta_solution,
            'n_eff': n_eff_solution,
            'V': self.V,
            'wavelength': self.wavelength,
            'n1': self.n1,
            'n2': self.n2,
            'a': self.a,
        }
    
    @property
    def u(self):
        """Normalized core parameter (requires solve_dispersion first)."""
        if self._u is None:
            self.solve_dispersion()
        return self._u
    
    @property
    def w(self):
        """Normalized cladding parameter."""
        if self._w is None:
            self.solve_dispersion()
        return self._w
    
    @property
    def beta(self):
        """Propagation constant β (µm⁻¹)."""
        if self._beta is None:
            self.solve_dispersion()
        return self._beta
    
    @property
    def n_eff(self):
        """Effective refractive index n_eff = β/k₀."""
        if self._n_eff is None:
            self.solve_dispersion()
        return self._n_eff
    
    def mode_profile(self, r):
        """
        Compute the LP₀₁ mode field profile Ψ(r).
        
        The mode profile is normalized such that:
            ∫₀^∞ |Ψ(r)|² 2π r dr = 1
        
        Parameters
        ----------
        r : float or ndarray
            Radial coordinate(s) in micrometers
        
        Returns
        -------
        ndarray
            Normalized mode profile Ψ(r)
        """
        r = np.atleast_1d(r)
        psi = np.zeros_like(r, dtype=float)
        
        # Core region (r < a)
        core_mask = r < self.a
        r_core = r[core_mask]
        if len(r_core) > 0:
            xi_core = self.u * r_core / self.a
            psi[core_mask] = jv(0, xi_core)
        
        # Cladding region (r >= a)
        clad_mask = r >= self.a
        r_clad = r[clad_mask]
        if len(r_clad) > 0:
            xi_clad = self.w * r_clad / self.a
            psi[clad_mask] = jv(0, self.u) / kv(0, self.w) * kv(0, xi_clad)
        
        # Normalize
        psi_normalized = psi / self._normalization_factor()
        
        return psi_normalized
    
    def _normalization_factor(self):
        """
        Compute the normalization constant for mode profile.
        
        Returns ∫₀^∞ |Ψ(r)|² 2π r dr (before normalization)
        """
        # Numerical integration using composite approach
        
        # Core region [0, a]
        r_core = np.linspace(0, self.a * 0.999, 200)
        xi_core = self.u * r_core / self.a
        psi_core = jv(0, xi_core)
        integrand_core = psi_core**2 * 2 * np.pi * r_core
        integral_core = np.trapz(integrand_core, r_core)
        
        # Cladding region [a, ∞)
        # Use exponential transform: r = a + s with dr = ds
        # K₀ decays exponentially, so we can integrate over finite domain
        r_clad = np.linspace(self.a, self.a + 10 / self.w, 300)
        xi_clad = self.w * r_clad / self.a
        psi_clad = jv(0, self.u) / kv(0, self.w) * kv(0, xi_clad)
        integrand_clad = psi_clad**2 * 2 * np.pi * r_clad
        integral_clad = np.trapz(integrand_clad, r_clad)
        
        return np.sqrt(integral_core + integral_clad)
    
    def summary(self):
        """Print a summary of fiber parameters and solution."""
        print("\n" + "="*60)
        print("STEP-INDEX FIBER PARAMETERS")
        print("="*60)
        print(f"Core index n₁:        {self.n1:.6f}")
        print(f"Cladding index n₂:    {self.n2:.6f}")
        print(f"Core radius a:        {self.a:.4f} µm")
        print(f"Wavelength λ₀:        {self.wavelength:.4f} µm")
        print(f"Wave number k₀:       {self.k0:.6f} µm⁻¹")
        print(f"V-number:             {self.V:.6f}")
        print("-"*60)
        print(f"Propagation β:        {self.beta:.8f} µm⁻¹")
        print(f"Effective index n_eff: {self.n_eff:.8f}")
        print(f"Normalized u:         {self.u:.8f}")
        print(f"Normalized w:         {self.w:.8f}")
        print("="*60 + "\n")


# Default test case: SMF-28-like fiber at telecom wavelength
if __name__ == "__main__":
    # Standard single-mode fiber parameters
    n1 = 1.48
    n2 = 1.46
    a = 4.1  # Core radius in µm
    wavelength = 1.55  # Telecom C-band in µm
    
    # Create fiber and solve
    fiber = StepIndexFiber(n1, n2, a, wavelength)
    solution = fiber.solve_dispersion(verbose=True)
    fiber.summary()
    
    # Test mode profile at a few points
    print("Testing mode profile computation...")
    r_test = np.array([0, 2.0, 4.1, 5.0, 10.0])
    psi = fiber.mode_profile(r_test)
    
    print(f"\nMode profile at test points:")
    print(f"{'r (µm)':>10} | {'Ψ(r)':>15} | {'Region':>10}")
    print("-"*40)
    for ri, psi_i in zip(r_test, psi):
        region = "core" if ri < a else "cladding"
        print(f"{ri:10.2f} | {psi_i:15.8f} | {region:>10}")
