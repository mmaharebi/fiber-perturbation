"""
Data Generation Script: Produce comprehensive numerical results for all scenarios.

This script generates validation data for:
1. Dispersion curves (n_eff vs wavelength)
2. LP₀₁ mode profiles (normalized field)
3. Radius perturbation validation (theory vs exact)
4. Index perturbation sensitivity maps
5. Absorption/attenuation vs conductivity

All results saved to /results/ for plotting.
"""

import numpy as np
import json
import pickle
import sys
sys.path.insert(0, '/home/mahdi/Desktop/fiber-perturbation/src')

from fiber_modes import StepIndexFiber
from perturbation import PerturbationAnalysis


class DataGenerator:
    """Generate comprehensive numerical results."""
    
    def __init__(self, output_dir='results'):
        self.output_dir = output_dir
        self.data_cache = {}
    
    def generate_dispersion_curve(self):
        """Generate effective index vs wavelength (dispersion curve)."""
        print("\n" + "="*70)
        print("GENERATING: Dispersion Curve (n_eff vs wavelength)")
        print("="*70)
        
        # Only use wavelengths where the fiber guides LP01 (V > 2.405)
        # For a=4.1µm and Δn=0.02, this means λ < ~2.5µm
        wavelengths = np.linspace(1.2, 1.8, 50)  # Reduced range to ensure V > 2.405
        n_eff_values = []
        beta_values = []
        V_values = []
        
        valid_wavelengths = []
        
        for wl in wavelengths:
            fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=wl)
            
            # Skip if fiber doesn't guide LP01 mode
            if fiber.V < 2.405:
                print(f"  Skipping λ={wl:.4f} µm (V={fiber.V:.4f} < 2.405 - below cutoff)")
                continue
            
            try:
                sol = fiber.solve_dispersion()
                n_eff_values.append(sol['n_eff'])
                beta_values.append(sol['beta'])
                V_values.append(sol['V'])
                valid_wavelengths.append(wl)
            except Exception as e:
                print(f"  Failed at λ={wl:.4f} µm: {e}")
                continue
        
        data = {
            'wavelength': valid_wavelengths,
            'n_eff': n_eff_values,
            'beta': beta_values,
            'V': V_values,
        }
        
        # Save
        with open(f'{self.output_dir}/dispersion_curve.json', 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✓ Generated {len(valid_wavelengths)} points")
        if valid_wavelengths:
            print(f"  Wavelength range: {valid_wavelengths[0]:.2f}-{valid_wavelengths[-1]:.2f} µm")
            print(f"  n_eff range: {min(n_eff_values):.6f}-{max(n_eff_values):.6f}")
        
        self.data_cache['dispersion'] = data
        return data
    
    def generate_mode_profiles(self):
        """Generate LP₀₁ mode radial profile."""
        print("\n" + "="*70)
        print("GENERATING: LP₀₁ Mode Profiles")
        print("="*70)
        
        fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
        fiber.solve_dispersion()
        
        # Generate profile over core and cladding
        r_core = np.linspace(0, fiber.a * 0.999, 200)
        r_clad = np.linspace(fiber.a, fiber.a + 15/fiber.w, 300)
        r_all = np.concatenate([r_core, r_clad])
        
        psi = fiber.mode_profile(r_all)
        
        data = {
            'r': r_all.tolist(),
            'psi': psi.tolist(),
            'psi_squared': (psi**2).tolist(),
            'core_radius': fiber.a,
            'wavelength': fiber.wavelength,
        }
        
        with open(f'{self.output_dir}/mode_profile.json', 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✓ Generated profile over r ∈ [0, {r_all[-1]:.2f}] µm")
        print(f"  Core radius: {fiber.a:.2f} µm")
        print(f"  Peak value: {np.max(np.abs(psi)):.6f}")
        
        self.data_cache['mode_profile'] = data
        return data
    
    def generate_radius_validation(self):
        """Generate radius perturbation validation data."""
        print("\n" + "="*70)
        print("GENERATING: Radius Perturbation Validation")
        print("="*70)
        
        fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
        fiber.solve_dispersion()
        perturb = PerturbationAnalysis(fiber)
        
        # Parameter range
        delta_a_range = np.linspace(-0.1, 0.1, 41)
        delta_a_range = delta_a_range[delta_a_range != 0]
        
        results = {
            'delta_a': [],
            'delta_a_percent': [],
            'delta_n_eff_theory': [],
            'delta_n_eff_exact': [],
            'relative_error': [],
        }
        
        for da in delta_a_range:
            result = perturb.scenario_2_radius_variation(da)
            results['delta_a'].append(da)
            results['delta_a_percent'].append(result['delta_a_relative'])
            results['delta_n_eff_theory'].append(result['delta_n_eff_theory'])
            results['delta_n_eff_exact'].append(result['delta_n_eff_exact'])
            results['relative_error'].append(result['relative_error'] * 100)
        
        data = {
            'delta_a': results['delta_a'],
            'delta_a_percent': results['delta_a_percent'],
            'delta_n_eff_theory': results['delta_n_eff_theory'],
            'delta_n_eff_exact': results['delta_n_eff_exact'],
            'relative_error': results['relative_error'],
        }
        
        with open(f'{self.output_dir}/radius_validation.json', 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✓ Generated {len(delta_a_range)} perturbation points")
        print(f"  Δa range: {min(delta_a_range):.4f} to {max(delta_a_range):.4f} µm")
        print(f"  Max error: {max(results['relative_error']):.4f}%")
        print(f"  Mean error: {np.mean(results['relative_error']):.4f}%")
        
        self.data_cache['radius_validation'] = data
        return data
    
    def generate_index_sensitivity(self):
        """Generate index perturbation sensitivity map."""
        print("\n" + "="*70)
        print("GENERATING: Index Perturbation Sensitivity Map")
        print("="*70)
        
        fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
        fiber.solve_dispersion()
        perturb = PerturbationAnalysis(fiber)
        
        # 2D parameter space
        delta_n_range = np.linspace(0.0001, 0.005, 20)
        delta_width_range = np.linspace(0.2, 3.0, 20)
        
        # Create meshgrid
        delta_n_grid = []
        delta_width_grid = []
        delta_n_eff_grid = []
        
        for dn in delta_n_range:
            for dw in delta_width_range:
                result = perturb.scenario_1_localized_index(dn, dw)
                delta_n_grid.append(dn)
                delta_width_grid.append(dw)
                delta_n_eff_grid.append(result['delta_n_eff'])
        
        data = {
            'delta_n': delta_n_grid,
            'delta_width': delta_width_grid,
            'delta_n_eff': delta_n_eff_grid,
            'delta_n_range': delta_n_range.tolist(),
            'delta_width_range': delta_width_range.tolist(),
        }
        
        with open(f'{self.output_dir}/index_sensitivity.json', 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✓ Generated {len(delta_n_grid)} sensitivity points")
        print(f"  Δn range: {delta_n_range[0]:.6f} to {delta_n_range[-1]:.6f}")
        print(f"  δ range: {delta_width_range[0]:.2f} to {delta_width_range[-1]:.2f} µm")
        
        self.data_cache['index_sensitivity'] = data
        return data
    
    def generate_absorption_data(self):
        """Generate absorption/attenuation vs conductivity."""
        print("\n" + "="*70)
        print("GENERATING: Absorption - Attenuation vs Conductivity")
        print("="*70)
        
        fiber = StepIndexFiber(n1=1.48, n2=1.46, a=4.1, wavelength=1.55)
        fiber.solve_dispersion()
        perturb = PerturbationAnalysis(fiber)
        
        # Conductivity range (log scale)
        sigma_range = np.logspace(-4, 2, 50)
        
        alpha_values = []
        alpha_dB_values = []
        lambda_att_values = []
        
        for sigma in sigma_range:
            result = perturb.scenario_3_weak_absorption(sigma)
            alpha = result['attenuation_alpha']
            alpha_dB = result['attenuation_dB_per_m']
            
            alpha_values.append(alpha)
            alpha_dB_values.append(alpha_dB)
            
            # Attenuation length (1/alpha in cm)
            if alpha > 1e-10:
                lambda_att = 1 / alpha / 1e4
            else:
                lambda_att = np.inf
            lambda_att_values.append(lambda_att)
        
        data = {
            'sigma': sigma_range.tolist(),
            'alpha': alpha_values,
            'alpha_dB_per_m': alpha_dB_values,
            'lambda_att': lambda_att_values,
        }
        
        with open(f'{self.output_dir}/absorption_data.json', 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✓ Generated {len(sigma_range)} conductivity points")
        print(f"  σ range: {sigma_range[0]:.2e} to {sigma_range[-1]:.2e} S/m")
        print(f"  α range: {min(alpha_values):.2e} to {max(alpha_values):.2e} 1/µm")
        
        self.data_cache['absorption'] = data
        return data
    
    def generate_fiber_schematic_data(self):
        """Generate data for step-index fiber schematic."""
        print("\n" + "="*70)
        print("GENERATING: Fiber Schematic Data")
        print("="*70)
        
        a = 4.1  # Core radius
        n1 = 1.48
        n2 = 1.46
        
        # Radial profile data
        r_plot = np.array([0, a, a, 10])
        n_plot = np.array([n1, n1, n2, n2])
        
        data = {
            'core_radius': a,
            'cladding_outer': 10.0,
            'n1': n1,
            'n2': n2,
            'r_profile': r_plot.tolist(),
            'n_profile': n_plot.tolist(),
        }
        
        with open(f'{self.output_dir}/fiber_schematic.json', 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✓ Schematic data created")
        print(f"  Core radius: {a:.2f} µm")
        print(f"  n1 = {n1}, n2 = {n2}")
        
        self.data_cache['schematic'] = data
        return data
    
    def generate_all(self):
        """Generate all data sets."""
        print("\n" + "█"*70)
        print("DATA GENERATION PIPELINE")
        print("█"*70)
        
        self.generate_dispersion_curve()
        self.generate_mode_profiles()
        self.generate_radius_validation()
        self.generate_index_sensitivity()
        self.generate_absorption_data()
        self.generate_fiber_schematic_data()
        
        print("\n" + "█"*70)
        print("DATA GENERATION COMPLETE ✓")
        print("█"*70)
        print(f"\nAll data saved to: {self.output_dir}/")
        print("Files created:")
        print("  - dispersion_curve.json")
        print("  - mode_profile.json")
        print("  - radius_validation.json")
        print("  - index_sensitivity.json")
        print("  - absorption_data.json")
        print("  - fiber_schematic.json")


if __name__ == "__main__":
    generator = DataGenerator(output_dir='/home/mahdi/Desktop/fiber-perturbation/results')
    generator.generate_all()
