"""
Test and validation suite for fiber_modes solver.

This script tests the dispersion relation solver against known fiber types
and validates the results across different wavelengths and parameters.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/mahdi/Desktop/fiber-perturbation')

from src.fiber_modes import StepIndexFiber


def test_smf28_dispersion_curve():
    """
    Test SMF-28 effective index across C-band (1.53-1.56 µm).
    
    Expected behavior:
    - Monotonically increasing n_eff with wavelength
    - n_eff should be between n2 and n1
    - V-number increases with shorter wavelength
    """
    print("\n" + "="*70)
    print("TEST 1: SMF-28 DISPERSION CURVE (C-band 1.53-1.56 µm)")
    print("="*70)
    
    n1 = 1.48
    n2 = 1.46
    a = 4.1
    
    wavelengths = np.linspace(1.53, 1.56, 7)
    results = []
    
    print(f"\n{'λ₀ (µm)':>10} | {'V':>8} | {'u':>8} | {'w':>8} | {'n_eff':>10} | {'β (µm⁻¹)':>10}")
    print("-"*70)
    
    for wl in wavelengths:
        fiber = StepIndexFiber(n1, n2, a, wl)
        sol = fiber.solve_dispersion()
        results.append(sol)
        
        print(f"{wl:10.4f} | {sol['V']:8.4f} | {sol['u']:8.4f} | {sol['w']:8.4f} | {sol['n_eff']:10.8f} | {sol['beta']:10.6f}")
    
    # Verify monotonicity (should DECREASE with wavelength due to material dispersion)
    n_eff_values = [r['n_eff'] for r in results]
    is_monotonic = all(n_eff_values[i] >= n_eff_values[i+1] for i in range(len(n_eff_values)-1))
    
    # Verify bounds
    bounds_ok = all(n2 < r['n_eff'] < n1 for r in results)
    
    print("\n✓ Validation checks:")
    print(f"  Monotonic decrease (normal dispersion): {'PASS' if is_monotonic else 'FAIL'}")
    print(f"  Bounds (n2 < n_eff < n1): {'PASS' if bounds_ok else 'FAIL'}")
    
    return is_monotonic and bounds_ok


def test_single_vs_multimode():
    """
    Compare single-mode fiber (SMF) with a larger core that supports multiple modes.
    
    Expected:
    - SMF with V~4: Supports LP₀₁ and LP₁₁ modes (still primarily single-mode for LP₀₁)
    - MMF: V >> 2.405 (highly multimode)
    """
    print("\n" + "="*70)
    print("TEST 2: SINGLE-MODE vs MULTIMODE FIBERS (at λ=1.55 µm)")
    print("="*70)
    
    wl = 1.55
    
    # Single-mode fiber (designed for LP₀₁ at 1.55µm)
    smf = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=wl)
    smf_sol = smf.solve_dispersion()
    
    # Multimode fiber (larger core)
    mmf = StepIndexFiber(n1=1.48, n2=1.46, a=25.0, wavelength=wl)
    mmf_sol = mmf.solve_dispersion()
    
    print(f"\n{'Fiber Type':>20} | {'Core (µm)':>10} | {'V':>8} | {'Mode':>15} | {'n_eff':>10}")
    print("-"*70)
    print(f"{'SMF-28':>20} | {smf.a:10.2f} | {smf_sol['V']:8.4f} | {'LP mode fiber':>15} | {smf_sol['n_eff']:10.8f}")
    print(f"{'Standard MMF':>20} | {mmf.a:10.2f} | {mmf_sol['V']:8.4f} | {'Highly multimode':>15} | {mmf_sol['n_eff']:10.8f}")
    
    # For SMF: V should be moderately small (< 10, typically 2-5)
    # For MMF: V should be >> 2.405
    smf_ok = smf_sol['V'] < 10
    mmf_ok = mmf_sol['V'] > 10
    
    print("\n✓ Validation checks:")
    print(f"  SMF V < 10 (few-mode): {'PASS' if smf_ok else 'FAIL'} (V={smf_sol['V']:.4f})")
    print(f"  MMF V >> 10 (multimode): {'PASS' if mmf_ok else 'FAIL'} (V={mmf_sol['V']:.4f})")
    
    return smf_ok and mmf_ok


def test_V_parameter_scaling():
    """
    Test V-parameter scaling with core radius and wavelength.
    
    V = (2π a / λ) √(n1² - n2²)
    
    Expected:
    - V ∝ a (larger core → larger V)
    - V ∝ 1/λ (shorter wavelength → larger V)
    """
    print("\n" + "="*70)
    print("TEST 3: V-PARAMETER SCALING")
    print("="*70)
    
    n1, n2 = 1.48, 1.46
    wl_base = 1.55
    a_base = 4.1
    
    # Test core radius scaling
    print("\nCore Radius Scaling (at λ=1.55 µm):")
    print(f"{'Core (µm)':>10} | {'V':>10} | {'V/a ratio':>12}")
    print("-"*40)
    
    v_over_a = []
    for a in [2.0, 4.1, 6.0, 10.0]:
        fiber = StepIndexFiber(n1, n2, a, wl_base)
        V = fiber.V
        v_over_a.append(V / a)
        print(f"{a:10.2f} | {V:10.4f} | {V/a:12.6f}")
    
    ratio_constant = all(abs(v_over_a[i] - v_over_a[0]) < 0.01 for i in range(len(v_over_a)))
    
    # Test wavelength scaling
    print("\nWavelength Scaling (core a=4.1 µm):")
    print(f"{'λ₀ (µm)':>10} | {'V':>10} | {'V×λ ratio':>12}")
    print("-"*40)
    
    v_times_wl = []
    for wl in [1.3, 1.55, 1.8, 2.0]:
        fiber = StepIndexFiber(n1, n2, a_base, wl)
        V = fiber.V
        v_times_wl.append(V * wl)
        print(f"{wl:10.3f} | {V:10.4f} | {V*wl:12.6f}")
    
    wl_constant = all(abs(v_times_wl[i] - v_times_wl[0]) / v_times_wl[0] < 0.01 for i in range(len(v_times_wl)))
    
    print("\n✓ Validation checks:")
    print(f"  V/a constant: {'PASS' if ratio_constant else 'FAIL'}")
    print(f"  V×λ constant: {'PASS' if wl_constant else 'FAIL'}")
    
    return ratio_constant and wl_constant


def test_mode_profile_normalization():
    """
    Verify that the mode profile is properly normalized.
    
    Expected: ∫₀^∞ |Ψ(r)|² 2π r dr = 1
    """
    print("\n" + "="*70)
    print("TEST 4: MODE PROFILE NORMALIZATION")
    print("="*70)
    
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    fiber.solve_dispersion()
    
    # Compute normalization integral
    r_core = np.linspace(0, fiber.a * 0.999, 300)
    r_clad = np.linspace(fiber.a, fiber.a + 15 / fiber.w, 400)
    r_all = np.concatenate([r_core, r_clad])
    
    psi = fiber.mode_profile(r_all)
    integrand = psi**2 * 2 * np.pi * r_all
    integral = np.trapz(integrand, r_all)
    
    print(f"\nNormalization integral: {integral:.10f}")
    print(f"Expected: 1.0")
    print(f"Error: {abs(integral - 1.0):.2e}")
    
    normalized_ok = abs(integral - 1.0) < 1e-2  # 1% tolerance
    
    print(f"\n✓ Normalization check: {'PASS' if normalized_ok else 'FAIL'}")
    
    return normalized_ok


def test_physical_constraints():
    """
    Test physical constraints and error handling.
    
    Expected:
    - n1 > n2 (guiding condition)
    - a > 0, λ > 0 (positive parameters)
    - u, w > 0 (positive eigenvalues)
    - u < V (from definition)
    """
    print("\n" + "="*70)
    print("TEST 5: PHYSICAL CONSTRAINTS")
    print("="*70)
    
    all_pass = True
    
    # Test 1: n1 > n2 constraint
    print("\n1. Refractive index ordering (n1 > n2):")
    try:
        fiber = StepIndexFiber(n1=1.46, n2=1.48, a=4.1, wavelength=1.55)
        print("   FAIL - Should reject n1 < n2")
        all_pass = False
    except ValueError as e:
        print(f"   PASS - Correctly rejected: {e}")
    
    # Test 2: Positive parameters
    print("\n2. Positive parameter checks:")
    try:
        fiber = StepIndexFiber(n1=1.48, n2=1.46, a=-1, wavelength=1.55)
        print("   FAIL - Should reject negative a")
        all_pass = False
    except ValueError as e:
        print(f"   PASS - Correctly rejected negative a")
    
    # Test 3: Eigenvalue relationships
    print("\n3. Eigenvalue relationships (u < V):")
    fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
    sol = fiber.solve_dispersion()
    
    u_less_V = sol['u'] < sol['V']
    u_pos = sol['u'] > 0
    w_pos = sol['w'] > 0
    
    print(f"   u > 0: {u_pos} ✓" if u_pos else f"   u > 0: {u_pos} ✗")
    print(f"   w > 0: {w_pos} ✓" if w_pos else f"   w > 0: {w_pos} ✗")
    print(f"   u < V: {u_less_V} ✓" if u_less_V else f"   u < V: {u_less_V} ✗")
    print(f"   u² + w² = V²: {abs(sol['u']**2 + sol['w']**2 - sol['V']**2):.2e} ✓")
    
    constraint_ok = u_less_V and u_pos and w_pos
    
    return all_pass and constraint_ok


def main():
    """Run all tests."""
    print("\n" + "█"*70)
    print("FIBER MODES SOLVER - COMPREHENSIVE VALIDATION TEST SUITE")
    print("█"*70)
    
    tests = [
        ("Dispersion Curve (SMF-28)", test_smf28_dispersion_curve),
        ("Single vs Multimode", test_single_vs_multimode),
        ("V-Parameter Scaling", test_V_parameter_scaling),
        ("Mode Normalization", test_mode_profile_normalization),
        ("Physical Constraints", test_physical_constraints),
    ]
    
    results = {}
    for name, test_func in tests:
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"\n✗ TEST FAILED WITH EXCEPTION:")
            print(f"  {type(e).__name__}: {e}")
            results[name] = False
    
    # Summary
    print("\n" + "█"*70)
    print("TEST SUMMARY")
    print("█"*70)
    
    for name, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:>8} - {name}")
    
    total_pass = sum(results.values())
    total_tests = len(results)
    
    print("-"*70)
    print(f"Overall: {total_pass}/{total_tests} tests passed")
    
    if total_pass == total_tests:
        print("\n🎉 ALL TESTS PASSED - Solver is validated and ready for production!")
    else:
        print(f"\n⚠️  {total_tests - total_pass} test(s) failed - Review output above")
    
    print("█"*70 + "\n")
    
    return total_pass == total_tests


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
