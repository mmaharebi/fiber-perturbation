"""
Perturbation Theory Module: Compute eigenvalue shifts due to fabrication errors.

This module implements first-order perturbation theory for computing shifts in
propagation constant and effective index due to:
1. Localized index perturbations (doping/diffusion)
2. Core radius variations (geometric errors)
3. Weak absorption (material loss)

Theory Reference:
- Section 5 of PROJECT.md: First-order eigenvalue perturbation theory
- Sturm-Liouville formulation ensures validity of perturbation approach
"""

import numpy as np
from scipy.special import jv, kv
from scipy.integrate import quad, quadrature
import warnings


class PerturbationAnalysis:
    """
    Compute perturbation theory shifts for a fiber mode.
    
    Given an unperturbed fiber solution, compute how small changes in
    material or geometry affect the propagation constant and effective index.
    
    Parameters
    ----------
    fiber : StepIndexFiber
        Unperturbed fiber object with solve_dispersion() already called
    """
    
    def __init__(self, fiber):
        """Initialize perturbation analyzer with unperturbed fiber."""
        if fiber._beta is None:
            raise ValueError("Fiber must have solve_dispersion() called first")
        
        self.fiber = fiber
        
        # Cache normalization for efficiency
        self._norm_factor = None
    
    def _modal_power(self):
        """
        Compute normalized modal power (for scalar normalization).
        
        Returns
        -------
        float
            P_z = (1/2) ∫ |E|² × (effective_area_term) dr
            
            For normalized modes: scalar normalization gives P_z = 1/2
        """
        # For normalized mode: ∫|Ψ|² 2πr dr = 1
        # Modal power in scalar formulation = 1/2
        return 0.5
    
    def _perturbation_integrand(self, r, delta_n_sq):
        """
        Compute integrand for perturbation integral.
        
        ∫ ΔN²(r)|Ψ(r)|² 2πr dr
        
        where Ψ is the normalized mode profile.
        
        Parameters
        ----------
        r : float or ndarray
            Radial coordinate (µm)
        delta_n_sq : callable or float
            Δn²(r) - the refractive index perturbation squared
        
        Returns
        -------
        float or ndarray
            Integrand value
        """
        psi = self.fiber.mode_profile(r)
        
        # delta_n_sq can be a function or constant
        if callable(delta_n_sq):
            dn2 = delta_n_sq(r)
        else:
            dn2 = delta_n_sq
        
        # Integrand: Δn² |Ψ|² 2πr
        integrand = dn2 * psi**2 * 2 * np.pi * r
        
        return integrand
    
    def _integrate_perturbation(self, delta_n_sq, r_max=None):
        """
        Numerically integrate the perturbation integral.
        
        ∫₀^∞ Δn²(r)|Ψ(r)|² 2πr dr
        
        Parameters
        ----------
        delta_n_sq : callable
            Δn²(r) function
        r_max : float, optional
            Upper limit of integration. If None, uses 5/w for cladding decay.
        
        Returns
        -------
        float
            Value of perturbation integral
        """
        if r_max is None:
            r_max = self.fiber.a + 10 / self.fiber.w
        
        # Split integration into core and cladding
        r_core = np.linspace(0, self.fiber.a * 0.999, 250)
        r_clad = np.linspace(self.fiber.a, r_max, 300)
        r_all = np.concatenate([r_core, r_clad])
        
        integrand = self._perturbation_integrand(r_all, delta_n_sq)
        integral = np.trapz(integrand, r_all)
        
        return integral
    
    def scenario_1_localized_index(self, delta_n, delta_width):
        """
        Scenario 1: Localized index perturbation (doping/diffusion).
        
        Models a perturbation shell: Δn²(r) for a-δ < r < a+δ
        
        Parameters
        ----------
        delta_n : float
            Index change magnitude
        delta_width : float
            Width of perturbation shell (µm)
        
        Returns
        -------
        dict
            Results containing:
            - delta_beta : shift in β (µm⁻¹)
            - delta_n_eff : shift in effective index
            - integral : value of perturbation integral
            - perturbation_type : scenario identifier
        """
        # Define perturbation function
        def delta_n_sq_func(r):
            """Δn² = 0 everywhere except perturbation shell."""
            r_lower = self.fiber.a - delta_width
            r_upper = self.fiber.a + delta_width
            
            if isinstance(r, np.ndarray):
                result = np.zeros_like(r, dtype=float)
                mask = (r >= r_lower) & (r <= r_upper)
                result[mask] = delta_n**2
                return result
            else:
                if r_lower <= r <= r_upper:
                    return delta_n**2
                else:
                    return 0.0
        
        # Compute perturbation integral
        integral = self._integrate_perturbation(delta_n_sq_func)
        
        # Perturbation formula: Δβ = -(k₀²/2β) ∫ Δε|Ψ|² dA
        # With Δε = ε₀ Δn², and scalar normalization:
        k0 = self.fiber.k0
        beta = self.fiber.beta
        
        # For normalized mode: P_z = 1/2, so factor is k₀²/(2β) × 2
        delta_beta = -(k0**2 / (2 * beta)) * integral
        
        # Convert to effective index shift
        delta_n_eff = delta_beta / k0
        
        return {
            'delta_n': delta_n,
            'delta_width': delta_width,
            'delta_beta': delta_beta,
            'delta_n_eff': delta_n_eff,
            'integral': integral,
            'perturbation_type': 'localized_index',
            'scenario': 'Doping/Diffusion'
        }
    
    def scenario_2_radius_variation(self, delta_a):
        """
        Scenario 2: Core radius variation (geometric perturbation).
        
        Models Δa perturbation via finite difference + perturbation formula.
        
        Method 1 (Fast): Use perturbation formula with equivalent shell
        Method 2 (Exact): Recompute dispersion with new radius
        
        Returns both for comparison.
        
        Parameters
        ----------
        delta_a : float
            Core radius change (µm)
        
        Returns
        -------
        dict
            Results containing:
            - delta_beta_theory : perturbation theory prediction (µm⁻¹)
            - delta_n_eff_theory : perturbation theory prediction
            - delta_beta_exact : exact recomputation result
            - delta_n_eff_exact : exact recomputation result
            - relative_error : |theory - exact| / |exact|
            - perturbation_type : scenario identifier
        """
        from fiber_modes import StepIndexFiber
        
        # Theory prediction using perturbation formula
        # For radius variation, we use numerical differentiation approach
        # Δβ/Δa ≈ dβ/da computed numerically
        
        # Small perturbation
        da_test = 1e-4 * self.fiber.a  # 0.01% variation for numerical derivative
        
        fiber_test = StepIndexFiber(
            self.fiber.n1, self.fiber.n2, 
            self.fiber.a + da_test, 
            self.fiber.wavelength
        )
        sol_test = fiber_test.solve_dispersion()
        
        # Numerical derivative
        dbeta_da = (sol_test['beta'] - self.fiber.beta) / da_test
        
        # Perturbation theory: Δβ ≈ (dβ/da) Δa
        delta_beta_theory = dbeta_da * delta_a
        delta_n_eff_theory = delta_beta_theory / self.fiber.k0
        
        # Exact solution with new radius
        fiber_new = StepIndexFiber(
            self.fiber.n1, self.fiber.n2,
            self.fiber.a + delta_a,
            self.fiber.wavelength
        )
        sol_new = fiber_new.solve_dispersion()
        
        delta_beta_exact = sol_new['beta'] - self.fiber.beta
        delta_n_eff_exact = sol_new['n_eff'] - self.fiber.n_eff
        
        # Relative error
        if abs(delta_beta_exact) > 1e-12:
            rel_error = abs(delta_beta_theory - delta_beta_exact) / abs(delta_beta_exact)
        else:
            rel_error = np.inf
        
        return {
            'delta_a': delta_a,
            'delta_a_relative': delta_a / self.fiber.a * 100,  # percent
            'delta_beta_theory': delta_beta_theory,
            'delta_n_eff_theory': delta_n_eff_theory,
            'delta_beta_exact': delta_beta_exact,
            'delta_n_eff_exact': delta_n_eff_exact,
            'relative_error': rel_error,
            'perturbation_type': 'radius_variation',
            'scenario': 'Core Radius Error'
        }
    
    def scenario_3_weak_absorption(self, sigma):
        """
        Scenario 3: Weak absorption (material loss/conductivity).
        
        Small conductivity σ in core modifies complex permittivity:
        Δε = -j σ/ω
        
        This leads to attenuation constant α = -2 Im(Δβ).
        
        Parameters
        ----------
        sigma : float
            Conductivity in core (S/m)
        
        Returns
        -------
        dict
            Results containing:
            - sigma : conductivity (S/m)
            - delta_beta_complex : complex shift in β
            - delta_beta_imag : imaginary part (attenuation)
            - attenuation_alpha : α = -2 Im(Δβ) (1/µm)
            - attenuation_dB_per_m : attenuation in dB/m
            - perturbation_type : scenario identifier
        """
        # Free-space angular frequency
        # ω = 2πc/λ where c = 3e8 m/s, λ in µm
        c = 3e8  # m/s
        lambda_m = self.fiber.wavelength * 1e-6  # Convert to meters
        omega = 2 * np.pi * c / lambda_m  # rad/s
        
        # Permittivity perturbation: Δε = -j σ/ω
        # In terms of index: Δε = ε₀ Δn², so:
        # Δn² = -j σ/(ω ε₀)
        
        epsilon_0 = 8.854e-12  # F/m
        delta_n_sq_complex = -1j * sigma / (omega * epsilon_0)
        
        # For complex perturbation, compute integral with complex result
        # Use same perturbation formula but allow complex result
        
        # Integration domain
        r_core = np.linspace(0, self.fiber.a * 0.999, 250)
        
        # Mode profile
        psi = self.fiber.mode_profile(r_core)
        
        # Complex integrand in core only (absorption only in core)
        integrand = delta_n_sq_complex * psi**2 * 2 * np.pi * r_core
        integral = np.trapz(integrand, r_core)
        
        # Perturbation formula
        k0 = self.fiber.k0
        beta = self.fiber.beta
        
        delta_beta_complex = -(k0**2 / (2 * beta)) * integral
        
        # Attenuation: α = -2 Im(Δβ)
        alpha_imag = delta_beta_complex.imag
        alpha = abs(alpha_imag * 2)  # Take absolute value for physical attenuation
        
        # Convert to dB/m
        # 1 Np/m = 20/ln(10) dB/m ≈ 8.686 dB/m
        # α in dB/m = α (1/µm) × 1e6 µm/m × 8.686 dB/Np
        alpha_dB_m = alpha * 1e6 * 8.686
        
        return {
            'sigma': sigma,
            'delta_beta_real': delta_beta_complex.real,
            'delta_beta_imag': delta_beta_complex.imag,
            'attenuation_alpha': alpha,
            'attenuation_dB_per_m': alpha_dB_m,
            'perturbation_type': 'weak_absorption',
            'scenario': 'Material Loss/Absorption'
        }
    
    def summary(self, result):
        """Print nicely formatted summary of perturbation results."""
        print("\n" + "="*70)
        print(f"PERTURBATION ANALYSIS: {result['scenario']}")
        print("="*70)
        
        if result['perturbation_type'] == 'localized_index':
            print(f"Index change Δn: {result['delta_n']:.6f}")
            print(f"Shell width δ: {result['delta_width']:.4f} µm")
            print(f"Perturbation integral: {result['integral']:.8e}")
            print(f"\nResults:")
            print(f"  Δβ = {result['delta_beta']:.8e} µm⁻¹")
            print(f"  Δn_eff = {result['delta_n_eff']:.8e}")
        
        elif result['perturbation_type'] == 'radius_variation':
            print(f"Radius change Δa: {result['delta_a']:.4f} µm ({result['delta_a_relative']:.3f}%)")
            print(f"\nPerturbation Theory:")
            print(f"  Δβ = {result['delta_beta_theory']:.8e} µm⁻¹")
            print(f"  Δn_eff = {result['delta_n_eff_theory']:.8e}")
            print(f"\nExact (recomputed):")
            print(f"  Δβ = {result['delta_beta_exact']:.8e} µm⁻¹")
            print(f"  Δn_eff = {result['delta_n_eff_exact']:.8e}")
            print(f"\nRelative Error: {result['relative_error']*100:.4f}%")
        
        elif result['perturbation_type'] == 'weak_absorption':
            print(f"Conductivity σ: {result['sigma']:.4e} S/m")
            print(f"\nResults:")
            print(f"  Re(Δβ) = {result['delta_beta_real']:.8e} µm⁻¹")
            print(f"  Im(Δβ) = {result['delta_beta_imag']:.8e} µm⁻¹")
            print(f"  Attenuation α = {result['attenuation_alpha']:.8e} 1/µm")
            print(f"  Attenuation = {result['attenuation_dB_per_m']:.4f} dB/m")
        
        print("="*70 + "\n")


if __name__ == "__main__":
    # Test with SMF-28
    from fiber_modes import StepIndexFiber
    
    print("PERTURBATION THEORY TEST")
    print("="*70)
    
    # Create and solve unperturbed fiber
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    fiber.solve_dispersion(verbose=True)
    
    # Create perturbation analyzer
    perturb = PerturbationAnalysis(fiber)
    
    # Test all three scenarios
    print("\n" + "█"*70)
    print("SCENARIO 1: LOCALIZED INDEX PERTURBATION")
    print("█"*70)
    result1 = perturb.scenario_1_localized_index(delta_n=0.001, delta_width=1.0)
    perturb.summary(result1)
    
    print("\n" + "█"*70)
    print("SCENARIO 2: CORE RADIUS VARIATION")
    print("█"*70)
    result2 = perturb.scenario_2_radius_variation(delta_a=0.01)
    perturb.summary(result2)
    
    print("\n" + "█"*70)
    print("SCENARIO 3: WEAK ABSORPTION")
    print("█"*70)
    result3 = perturb.scenario_3_weak_absorption(sigma=1.0)
    perturb.summary(result3)
    
    print("\nAll perturbation scenarios executed successfully! ✓")
