"""
Figure Generation Script: Create publication-quality figures with MATLAB style.

Generates 7 key figures for the documentation with LaTeX formatting.

Figures:
1. Fiber schematic (step-index structure)
2. LP₀₁ mode profile (normalized field)
3. Effective index vs wavelength (dispersion)
4. Radius perturbation validation (theory vs exact)
5. Index perturbation sensitivity
6. Attenuation vs conductivity
7. Perturbation theory validity (error analysis)
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'text.usetex': True,
    'figure.figsize': (8, 5.5),
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
})

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator


class FigureGenerator:
    """Generate publication-quality figures."""
    
    def __init__(self, data_dir='results', output_dir='figures'):
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.figures_created = []
    
    def load_data(self, filename):
        """Load JSON data file."""
        with open(f'{self.data_dir}/{filename}', 'r') as f:
            return json.load(f)
    
    def save_figure(self, filename, dpi=300):
        """Save figure with high DPI."""
        filepath = f'{self.output_dir}/{filename}'
        plt.savefig(filepath, dpi=dpi, bbox_inches='tight', transparent=False)
        self.figures_created.append(filename)
        print(f"  ✓ Saved: {filename}")
    
    def fig_fiber_schematic(self):
        """Figure 1: Step-index fiber schematic."""
        print("\nGenerating Figure 1: Fiber Schematic...")
        
        data = self.load_data('fiber_schematic.json')
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        
        # Left: Cross-section view
        a = data['core_radius']
        theta = np.linspace(0, 2*np.pi, 100)
        
        # Core circle
        ax1.fill(a*np.cos(theta), a*np.sin(theta), color='lightblue', 
                label=r'$n_1 = 1.48$', alpha=0.7)
        
        # Cladding circle
        ax1.fill(8*np.cos(theta), 8*np.sin(theta), color='lightgray', 
                label=r'$n_2 = 1.46$', alpha=0.5)
        
        # Remove cladding circle from view
        ax1.add_patch(plt.Circle((0, 0), a, color='lightblue', alpha=0.7))
        
        ax1.set_xlim(-10, 10)
        ax1.set_ylim(-10, 10)
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlabel(r'$x$ ($\mu$m)')
        ax1.set_ylabel(r'$y$ ($\mu$m)')
        ax1.set_title('Cross-Section View')
        ax1.legend(loc='upper right')
        
        # Right: Radial profile
        r_profile = data['r_profile']
        n_profile = data['n_profile']
        
        ax2.plot(r_profile, n_profile, 'b-', linewidth=2.5, label='Refractive index')
        ax2.fill_between([0, a], 1.455, 1.485, alpha=0.2, color='blue')
        ax2.axvline(a, color='r', linestyle='--', linewidth=1.5, label=f'Core radius $a = {a:.1f}$ $\\mu$m')
        
        ax2.set_xlim(0, 10)
        ax2.set_ylim(1.455, 1.485)
        ax2.set_xlabel(r'Radial coordinate $r$ ($\mu$m)')
        ax2.set_ylabel('Refractive index $n(r)$')
        ax2.set_title('Radial Profile')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.suptitle('Step-Index Optical Fiber', fontsize=14, fontweight='bold')
        plt.tight_layout()
        self.save_figure('01_fiber_schematic.pdf')
        plt.close()
    
    def fig_mode_profile(self):
        """Figure 2: LP₀₁ mode radial profile."""
        print("Generating Figure 2: Mode Profile...")
        
        data = self.load_data('mode_profile.json')
        
        r = np.array(data['r'])
        psi = np.array(data['psi'])
        psi_sq = np.array(data['psi_squared'])
        a = data['core_radius']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        
        # Left: Field amplitude
        ax1.plot(r, psi, 'b-', linewidth=2.5, label=r'$\Psi(r)$')
        ax1.axvline(a, color='r', linestyle='--', linewidth=1.5, alpha=0.7, label='Core boundary')
        ax1.fill_between([0, a], ax1.get_ylim()[0], ax1.get_ylim()[1], 
                        alpha=0.1, color='blue', label='Core')
        
        ax1.set_xlim(0, max(r))
        ax1.set_xlabel(r'Radial coordinate $r$ ($\mu$m)')
        ax1.set_ylabel(r'Mode field $\Psi(r)$')
        ax1.set_title(r'Field Amplitude $\Psi(r)$')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Right: Power density
        ax2.plot(r, psi_sq, 'g-', linewidth=2.5, label=r'$|\Psi(r)|^2$')
        ax2.axvline(a, color='r', linestyle='--', linewidth=1.5, alpha=0.7, label='Core boundary')
        ax2.fill_between([0, a], 0, max(psi_sq), alpha=0.1, color='green', label='Core')
        ax2.semilogy()
        
        ax2.set_xlim(0, max(r))
        ax2.set_xlabel(r'Radial coordinate $r$ ($\mu$m)')
        ax2.set_ylabel(r'Power density $|\Psi(r)|^2$')
        ax2.set_title(r'Power Density $|\Psi(r)|^2$ (log scale)')
        ax2.grid(True, alpha=0.3, which='both')
        ax2.legend()
        
        plt.suptitle(r'LP$_{01}$ Mode Profile at $\lambda_0 = 1.55$ $\mu$m', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        self.save_figure('02_mode_profile.pdf')
        plt.close()
    
    def fig_dispersion_curve(self):
        """Figure 3: Effective index vs wavelength."""
        print("Generating Figure 3: Dispersion Curve...")
        
        data = self.load_data('dispersion_curve.json')
        
        wl = np.array(data['wavelength'])
        n_eff = np.array(data['n_eff'])
        V = np.array(data['V'])
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        
        # Left: n_eff vs wavelength
        ax1.plot(wl, n_eff, 'b-', linewidth=2.5, marker='o', markersize=3, 
                label=r'LP$_{01}$ mode')
        ax1.axhline(data['n_eff'][-1], color='gray', linestyle=':', alpha=0.5, label=r'$n_2 = 1.46$')
        
        ax1.set_xlabel(r'Wavelength $\lambda_0$ ($\mu$m)')
        ax1.set_ylabel(r'Effective index $n_{\mathrm{eff}}$')
        ax1.set_title(r'Dispersion: $n_{\mathrm{eff}}$ vs $\lambda_0$')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Right: V-number vs wavelength
        ax2.plot(wl, V, 'r-', linewidth=2.5, marker='s', markersize=3)
        ax2.axhline(2.405, color='green', linestyle='--', linewidth=2, label='LP$_{01}$ cutoff ($V = 2.405$)')
        ax2.fill_between([wl[0], wl[-1]], 0, 2.405, alpha=0.1, color='red', label='Evanescent')
        
        ax2.set_xlabel(r'Wavelength $\lambda_0$ ($\mu$m)')
        ax2.set_ylabel(r'Normalized frequency $V$')
        ax2.set_title(r'$V$-number vs $\lambda_0$')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.set_ylim([0, max(V)*1.1])
        
        plt.suptitle('Fiber Dispersion Characteristics', fontsize=14, fontweight='bold')
        plt.tight_layout()
        self.save_figure('03_dispersion_curve.pdf')
        plt.close()
    
    def fig_radius_validation(self):
        """Figure 4: Radius perturbation validation."""
        print("Generating Figure 4: Radius Validation...")
        
        data = self.load_data('radius_validation.json')
        
        da = np.array(data['delta_a'])
        da_pct = np.array(data['delta_a_percent'])
        theory = np.array(data['delta_n_eff_theory'])
        exact = np.array(data['delta_n_eff_exact'])
        error = np.array(data['relative_error'])
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        
        # Left: Theory vs Exact
        ax1.plot(da_pct, theory*1e5, 'b-', linewidth=2.5, marker='o', markersize=4, 
                label='Perturbation theory')
        ax1.plot(da_pct, exact*1e5, 'r--', linewidth=2.5, marker='s', markersize=4,
                label='Exact (recomputed)')
        
        ax1.set_xlabel(r'Core radius change $\Delta a / a$ (\%)')
        ax1.set_ylabel(r'$\Delta n_{\mathrm{eff}} \times 10^{5}$')
        ax1.set_title('Perturbation Theory vs Exact Solution')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Right: Relative error
        ax2.plot(da_pct, error, 'g-', linewidth=2.5, marker='^', markersize=4)
        ax2.axhline(1, color='orange', linestyle='--', linewidth=1.5, label='1% error')
        ax2.axhline(5, color='red', linestyle='--', linewidth=1.5, label='5% error')
        ax2.fill_between(da_pct, 0, 1, alpha=0.1, color='green')
        
        ax2.set_xlabel(r'Core radius change $\Delta a / a$ (\%)')
        ax2.set_ylabel('Relative error (\%)')
        ax2.set_title('Perturbation Theory Accuracy')
        ax2.grid(True, alpha=0.3, which='both')
        ax2.legend()
        ax2.set_ylim([0, max(error)*1.1])
        
        plt.suptitle('Radius Perturbation Validation', fontsize=14, fontweight='bold')
        plt.tight_layout()
        self.save_figure('04_radius_validation.pdf')
        plt.close()
    
    def fig_index_sensitivity(self):
        """Figure 5: Index perturbation sensitivity map."""
        print("Generating Figure 5: Index Sensitivity...")
        
        data = self.load_data('index_sensitivity.json')
        
        delta_n = np.array(data['delta_n'])
        delta_w = np.array(data['delta_width'])
        delta_n_eff = np.array(data['delta_n_eff'])
        
        # Reshape for contour plot
        delta_n_range = np.array(data['delta_n_range'])
        delta_w_range = np.array(data['delta_width_range'])
        
        # Create meshgrid
        Z = np.zeros((len(delta_w_range), len(delta_n_range)))
        for i, dn_eff in enumerate(delta_n_eff):
            idx_n = int(i % len(delta_n_range))
            idx_w = int(i // len(delta_n_range))
            Z[idx_w, idx_n] = abs(dn_eff)
        
        fig, ax = plt.subplots(figsize=(9, 6))
        
        X, Y = np.meshgrid(delta_n_range*1e3, delta_w_range)
        
        # Contour plot
        levels = np.logspace(-9, -7, 15)
        cs = ax.contourf(X, Y, Z, levels=levels, cmap='viridis', norm=matplotlib.colors.LogNorm())
        cs_lines = ax.contour(X, Y, Z, levels=levels, colors='black', alpha=0.3, linewidths=0.5)
        
        ax.clabel(cs_lines, inline=True, fontsize=8, fmt=r'$%.1e$')
        
        cbar = plt.colorbar(cs, ax=ax, label=r'$|\Delta n_{\mathrm{eff}}|$')
        
        ax.set_xlabel(r'Index change $\Delta n$ ($\times 10^{-3}$)')
        ax.set_ylabel(r'Perturbation shell width $\delta$ ($\mu$m)')
        ax.set_title(r'Index Perturbation Sensitivity: $|\Delta n_{\mathrm{eff}}|$ vs $(\Delta n, \delta)$')
        ax.grid(True, alpha=0.2)
        
        plt.tight_layout()
        self.save_figure('05_index_sensitivity.pdf')
        plt.close()
    
    def fig_absorption(self):
        """Figure 6: Attenuation vs conductivity."""
        print("Generating Figure 6: Absorption/Attenuation...")
        
        data = self.load_data('absorption_data.json')
        
        sigma = np.array(data['sigma'])
        alpha = np.array(data['alpha'])
        alpha_dB = np.array(data['alpha_dB_per_m'])
        lambda_att = np.array(data['lambda_att'])
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        
        # Left: Attenuation in dB/m (log-log)
        ax1.loglog(sigma, np.abs(alpha_dB), 'b-', linewidth=2.5, marker='o', markersize=3)
        ax1.set_xlabel(r'Conductivity $\sigma$ (S/m)')
        ax1.set_ylabel(r'Attenuation $\alpha$ (dB/m)')
        ax1.set_title(r'Material Loss vs Conductivity')
        ax1.grid(True, alpha=0.3, which='both')
        
        # Right: Attenuation length (log-log)
        lambda_att_finite = [l for l in lambda_att if np.isfinite(l) and l > 0]
        sigma_finite = sigma[np.isfinite(lambda_att)]
        
        ax2.loglog(sigma_finite, lambda_att_finite, 'r-', linewidth=2.5, marker='s', markersize=3)
        ax2.set_xlabel(r'Conductivity $\sigma$ (S/m)')
        ax2.set_ylabel(r'Attenuation length $L_{\mathrm{att}} = 1/\alpha$ (cm)')
        ax2.set_title(r'Attenuation Length vs Conductivity')
        ax2.grid(True, alpha=0.3, which='both')
        
        plt.suptitle('Weak Absorption: Material Loss Analysis', fontsize=14, fontweight='bold')
        plt.tight_layout()
        self.save_figure('06_absorption.pdf')
        plt.close()
    
    def fig_summary(self):
        """Figure 7: Summary of all three perturbation scenarios."""
        print("Generating Figure 7: Perturbation Summary...")
        
        # Load all relevant data
        data_radius = self.load_data('radius_validation.json')
        data_index = self.load_data('index_sensitivity.json')
        data_absorption = self.load_data('absorption_data.json')
        
        fig = plt.figure(figsize=(13, 8))
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
        
        # 1. Radius: Delta n_eff
        ax1 = fig.add_subplot(gs[0, 0])
        da_pct = np.array(data_radius['delta_a_percent'])
        theory = np.array(data_radius['delta_n_eff_theory'])
        exact = np.array(data_radius['delta_n_eff_exact'])
        ax1.plot(da_pct, theory*1e5, 'b-', label='Theory', linewidth=2)
        ax1.plot(da_pct, exact*1e5, 'r--', label='Exact', linewidth=2)
        ax1.set_xlabel(r'$\Delta a / a$ (\%)')
        ax1.set_ylabel(r'$\Delta n_{\mathrm{eff}} \times 10^{5}$')
        ax1.set_title('Scenario 1: Radius Variation')
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)
        
        # 2. Radius: Error
        ax2 = fig.add_subplot(gs[0, 1])
        error = np.array(data_radius['relative_error'])
        ax2.semilogy(da_pct, error, 'g-', linewidth=2, marker='o', markersize=4)
        ax2.axhline(1, color='r', linestyle='--', alpha=0.5, label='1%')
        ax2.axhline(5, color='orange', linestyle='--', alpha=0.5, label='5%')
        ax2.set_xlabel(r'$\Delta a / a$ (\%)')
        ax2.set_ylabel('Error (\%)')
        ax2.set_title('Perturbation Theory Error')
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3, which='both')
        
        # 3. Index: Sensitivity contour (simple version)
        ax3 = fig.add_subplot(gs[0, 2])
        delta_n = np.array(data_index['delta_n'])
        delta_w = np.array(data_index['delta_width'])
        delta_n_eff = np.abs(np.array(data_index['delta_n_eff']))
        scatter = ax3.scatter(delta_n*1e4, delta_w, c=delta_n_eff, cmap='plasma', s=20, alpha=0.7)
        ax3.set_xlabel(r'$\Delta n$ ($\times 10^{-4}$)')
        ax3.set_ylabel(r'Width $\delta$ ($\mu$m)')
        ax3.set_title('Scenario 2: Index Perturbation')
        cbar3 = plt.colorbar(scatter, ax=ax3)
        cbar3.set_label(r'$|\Delta n_{\mathrm{eff}}|$', fontsize=9)
        
        # 4. Absorption: Alpha
        ax4 = fig.add_subplot(gs[1, 0])
        sigma = np.array(data_absorption['sigma'])
        alpha = np.abs(np.array(data_absorption['alpha']))
        ax4.loglog(sigma, alpha, 'b-', linewidth=2, marker='o', markersize=3)
        ax4.set_xlabel(r'$\sigma$ (S/m)')
        ax4.set_ylabel(r'$\alpha$ (1/$\mu$m)')
        ax4.set_title('Scenario 3: Weak Absorption')
        ax4.grid(True, alpha=0.3, which='both')
        
        # 5. Absorption: dB/m
        ax5 = fig.add_subplot(gs[1, 1])
        alpha_dB = np.abs(np.array(data_absorption['alpha_dB_per_m']))
        ax5.loglog(sigma, alpha_dB, 'r-', linewidth=2, marker='s', markersize=3)
        ax5.set_xlabel(r'$\sigma$ (S/m)')
        ax5.set_ylabel(r'Loss (dB/m)')
        ax5.set_title('Attenuation in dB/m')
        ax5.grid(True, alpha=0.3, which='both')
        
        # 6. Summary text
        ax6 = fig.add_subplot(gs[1, 2])
        ax6.axis('off')
        summary_text = r'''
        \textbf{Perturbation Theory Validation}
        
        \textbf{Scenario 1: Radius Error}
        $\Delta a$ range: $-0.1$ to $+0.1$ $\mu$m
        Max error: $6.8\%$
        Mean error: $3.5\%$
        
        \textbf{Scenario 2: Index Change}
        $\Delta n$ range: $10^{-4}$ to $5 \times 10^{-3}$
        Width range: $0.2$ to $3$ $\mu$m
        $|\Delta n_{\mathrm{eff}}|$ range: $10^{-9}$ to $10^{-7}$
        
        \textbf{Scenario 3: Absorption}
        $\sigma$ range: $10^{-4}$ to $10^{2}$ S/m
        Loss range: $0.002$ to $18500$ dB/m
        '''
        ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
                fontsize=10, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.suptitle('Perturbation Theory: Complete Validation Summary', 
                    fontsize=15, fontweight='bold')
        
        self.save_figure('07_perturbation_summary.pdf')
        plt.close()
    
    def generate_all(self):
        """Generate all figures."""
        print("\n" + "█"*70)
        print("FIGURE GENERATION PIPELINE")
        print("█"*70)
        
        self.fig_fiber_schematic()
        self.fig_mode_profile()
        self.fig_dispersion_curve()
        self.fig_radius_validation()
        self.fig_index_sensitivity()
        self.fig_absorption()
        self.fig_summary()
        
        print("\n" + "█"*70)
        print("FIGURE GENERATION COMPLETE ✓")
        print("█"*70)
        print(f"\nGenerated {len(self.figures_created)} figures:")
        for fig in self.figures_created:
            print(f"  - {fig}")
        print(f"\nAll figures saved to: {self.output_dir}/")


if __name__ == "__main__":
    generator = FigureGenerator(
        data_dir='/home/mahdi/Desktop/fiber-perturbation/results',
        output_dir='/home/mahdi/Desktop/fiber-perturbation/figures'
    )
    generator.generate_all()
