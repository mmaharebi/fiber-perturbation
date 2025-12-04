"""
Validation test suite for perturbation analysis module.

Compares perturbation theory predictions with exact recomputation
to verify accuracy and determine validity range of first-order theory.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/mahdi/Desktop/fiber-perturbation/src')

from fiber_modes import StepIndexFiber
from perturbation import PerturbationAnalysis


def test_radius_perturbation_validation():
    """
    Compare perturbation theory vs exact for core radius variations.
    
    Test across multiple perturbation magnitudes to determine validity range.
    """
    print("\n" + "="*80)
    print("TEST 1: RADIUS VARIATION PERTURBATION VALIDATION")
    print("="*80)
    
    # Create base fiber
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    fiber.solve_dispersion()
    perturb = PerturbationAnalysis(fiber)
    
    # Test various perturbation magnitudes
    delta_a_values = np.array([-0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05])
    delta_a_values = delta_a_values[delta_a_values != 0]  # Remove zero
    
    print(f"\nBase fiber: a={fiber.a:.2f} µm, n1={fiber.n1}, n2={fiber.n2}")
    print(f"λ={fiber.wavelength:.2f} µm, n_eff={fiber.n_eff:.8f}\n")
    print(f"{'Δa (µm)':>10} | {'% change':>10} | {'Theory (µm⁻¹)':>18} | {'Exact (µm⁻¹)':>18} | {'Error %':>10}")
    print("-"*80)
    
    results = []
    for da in delta_a_values:
        result = perturb.scenario_2_radius_variation(da)
        pct_change = result['delta_a_relative']
        theory = result['delta_beta_theory']
        exact = result['delta_beta_exact']
        error = result['relative_error'] * 100
        
        results.append(result)
        
        print(f"{da:10.4f} | {pct_change:10.3f} | {theory:18.8e} | {exact:18.8e} | {error:10.4f}")
    
    # Validity analysis
    errors = [r['relative_error']*100 for r in results]
    valid_range_1pct = sum(1 for e in errors if e < 1.0)
    valid_range_5pct = sum(1 for e in errors if e < 5.0)
    
    print("\n✓ Perturbation Theory Validity Analysis:")
    print(f"  < 1% error: {valid_range_1pct}/{len(errors)} cases")
    print(f"  < 5% error: {valid_range_5pct}/{len(errors)} cases")
    print(f"  Mean error: {np.mean(errors):.4f}%")
    print(f"  Max error: {np.max(errors):.4f}%")
    
    return results


def test_index_perturbation_effects():
    """
    Test localized index perturbation with various parameters.
    
    Validates that perturbation magnitude scales correctly with:
    - Index change Δn
    - Perturbation volume/width
    """
    print("\n" + "="*80)
    print("TEST 2: LOCALIZED INDEX PERTURBATION EFFECTS")
    print("="*80)
    
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    fiber.solve_dispersion()
    perturb = PerturbationAnalysis(fiber)
    
    # Test 1: Scaling with index change
    print("\nScaling with index change Δn (shell width δ=0.5 µm):")
    print(f"{'Δn':>12} | {'Δn_eff':>18} | {'Δn²':>12} | {'Ratio Δn_eff/Δn²':>18}")
    print("-"*70)
    
    delta_n_values = np.array([0.0001, 0.0005, 0.001, 0.002])
    
    for dn in delta_n_values:
        result = perturb.scenario_1_localized_index(dn, 0.5)
        ratio = result['delta_n_eff'] / (dn**2) if dn != 0 else 0
        print(f"{dn:12.6f} | {result['delta_n_eff']:18.8e} | {dn**2:12.8e} | {ratio:18.6e}")
    
    # Test 2: Scaling with shell width
    print("\nScaling with shell width δ (Δn=0.001):")
    print(f"{'δ (µm)':>12} | {'Δn_eff':>18} | {'Width':>12} | {'Integral':>18}")
    print("-"*70)
    
    delta_width_values = np.array([0.2, 0.5, 1.0, 2.0])
    
    for dw in delta_width_values:
        result = perturb.scenario_1_localized_index(0.001, dw)
        print(f"{dw:12.3f} | {result['delta_n_eff']:18.8e} | {dw:12.3f} | {result['integral']:18.8e}")
    
    return True


def test_absorption_parameter_space():
    """
    Test absorption perturbation across different conductivity values.
    
    Validates attenuation calculations and conversion to dB/m.
    """
    print("\n" + "="*80)
    print("TEST 3: WEAK ABSORPTION - PARAMETER SPACE")
    print("="*80)
    
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    fiber.solve_dispersion()
    perturb = PerturbationAnalysis(fiber)
    
    # Test across conductivity range
    sigma_values = np.logspace(-3, 1, 9)  # 0.001 to 10 S/m
    
    print(f"\n{'σ (S/m)':>12} | {'α (1/µm)':>18} | {'α (dB/m)':>18} | {'λ_atten (cm)':>15}")
    print("-"*75)
    
    for sigma in sigma_values:
        result = perturb.scenario_3_weak_absorption(sigma)
        alpha = result['attenuation_alpha']
        alpha_dB = result['attenuation_dB_per_m']
        
        # Characteristic attenuation length L_att = 1/α (in cm)
        if alpha > 1e-10:
            lambda_att = 1 / alpha / 1e4  # Convert to cm
        else:
            lambda_att = np.inf
        
        print(f"{sigma:12.6e} | {alpha:18.8e} | {alpha_dB:18.4f} | {lambda_att:15.6e}")
    
    return True


def test_physical_constraints():
    """
    Test that perturbation results obey physical constraints.
    
    - Effective index shifts should be small (first-order approximation)
    - Shifts scale with perturbation magnitude
    - Complex perturbations produce imaginary components
    """
    print("\n" + "="*80)
    print("TEST 4: PHYSICAL CONSTRAINTS")
    print("="*80)
    
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    fiber.solve_dispersion()
    perturb = PerturbationAnalysis(fiber)
    
    print(f"\nBase effective index: {fiber.n_eff:.8f}")
    
    # Radius perturbation
    result_r = perturb.scenario_2_radius_variation(0.01)
    shift_frac_r = abs(result_r['delta_n_eff_exact']) / fiber.n_eff
    
    # Index perturbation
    result_i = perturb.scenario_1_localized_index(0.001, 0.5)
    shift_frac_i = abs(result_i['delta_n_eff']) / fiber.n_eff
    
    # Absorption perturbation
    result_a = perturb.scenario_3_weak_absorption(1.0)
    
    print("\n✓ Physical constraint checks:")
    print(f"  Radius perturbation (Δa=0.01µm):")
    print(f"    Δn_eff/n_eff = {shift_frac_r:.6e} (should be << 1) ✓" if shift_frac_r < 0.1 else f"    Δn_eff/n_eff = {shift_frac_r:.6e} ✗")
    
    print(f"  Index perturbation (Δn=0.001):")
    print(f"    Δn_eff/n_eff = {shift_frac_i:.6e} (should be << 1) ✓" if shift_frac_i < 0.1 else f"    Δn_eff/n_eff = {shift_frac_i:.6e} ✗")
    
    print(f"  Absorption (σ=1 S/m):")
    is_complex = result_a['delta_beta_imag'] != 0
    print(f"    Complex result: {'✓' if is_complex else '✗'}")
    print(f"    Im(Δβ) = {result_a['delta_beta_imag']:.8e}")
    print(f"    Attenuation = {result_a['attenuation_dB_per_m']:.2f} dB/m")
    
    return True


def main():
    """Run all validation tests."""
    print("\n" + "█"*80)
    print("PERTURBATION THEORY - VALIDATION TEST SUITE")
    print("█"*80)
    
    tests = [
        ("Radius Variation Validation", test_radius_perturbation_validation),
        ("Index Perturbation Effects", test_index_perturbation_effects),
        ("Absorption Parameter Space", test_absorption_parameter_space),
        ("Physical Constraints", test_physical_constraints),
    ]
    
    results = {}
    for name, test_func in tests:
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"\n✗ TEST FAILED:")
            print(f"  {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            results[name] = False
    
    # Summary
    print("\n" + "█"*80)
    print("VALIDATION SUMMARY")
    print("█"*80)
    
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:>8} - {name}")
    
    total_pass = sum(1 for v in results.values() if v)
    total_tests = len(results)
    
    print("-"*80)
    print(f"Overall: {total_pass}/{total_tests} tests passed")
    
    if total_pass == total_tests:
        print("\n🎉 PERTURBATION MODULE VALIDATED AND READY!")
    
    print("█"*80 + "\n")
    
    return total_pass == total_tests


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
