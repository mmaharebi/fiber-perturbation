# Perturbation Theory for Optical Fibers: A Numerical Framework

> **Blog Post Draft for Personal Website**  
> **Author:** Mahdy Mahareb  
> **Date:** December 2025  
> **GitHub Repository:** [fiber-perturbation](https://github.com/mmaharebi/fiber-perturbation)  
> **Full Report (PDF):** [37-page technical document](https://github.com/mmaharebi/fiber-perturbation/blob/master/docs/main.pdf)

---

## Introduction

How do small fabrication errors affect optical fiber performance? This question is critical for fiber optic design, sensor applications, and integrated photonics. In this project, I developed a complete numerical framework to predict how tiny deviations in refractive index, core radius, or material properties shift the propagation constant and effective index of guided modes in step-index optical fibers.

The framework combines **analytical electromagnetic theory** with **first-order perturbation theory** and **numerical root-finding methods**, implemented entirely in Python with full validation and documentation.

---

## The Physical Problem

### Step-Index Optical Fiber

An optical fiber consists of a cylindrical core (refractive index n₁) surrounded by cladding (refractive index n₂ < n₁). Light propagates along the fiber axis by total internal reflection, confined to discrete **guided modes**.

![Fiber Geometry](figures/01_fiber_schematic.pdf)
*Figure 1: Step-index fiber cross-section and radial refractive index profile*

### Key Questions

1. How does the **propagation constant β** depend on wavelength, core radius, and index contrast?
2. What happens when fabrication introduces **small errors**:
   - Localized index changes (doping diffusion)
   - Core radius variations (geometric tolerances)
   - Weak absorption (material impurities)

---

## Mathematical Framework

### Dispersion Relation

For the fundamental LP₀₁ mode, the electric field satisfies boundary conditions that lead to the **transcendental eigenvalue equation**:

$$\frac{J_1(u)}{u J_0(u)} + \frac{K_1(w)}{w K_0(w)} = 0$$

where:
- $u = a\sqrt{k_0^2 n_1^2 - \beta^2}$ (core parameter)
- $w = a\sqrt{\beta^2 - k_0^2 n_2^2}$ (cladding decay parameter)
- $V = \sqrt{u^2 + w^2}$ (normalized frequency)
- J₀, J₁: Bessel functions, K₀, K₁: modified Bessel functions

This is solved numerically using **Brent's root-finding method** with automatic bracket detection.

### Perturbation Theory

When the fiber is perturbed (Δε(r) change in permittivity), **first-order perturbation theory** predicts:

$$\Delta\beta = \frac{\omega}{4P_z} \int_0^\infty \Delta\epsilon(r) |\Psi(r)|^2 2\pi r \, dr$$

where Ψ(r) is the unperturbed mode profile and Pz is the modal power.

---

## Implementation

### Core Modules

**`fiber_modes.py`** (380 lines)
- Solves LP₀₁ dispersion relation
- Computes normalized mode profiles
- Achieves 10⁻¹² convergence residual

**`perturbation.py`** (410 lines)
- Implements three perturbation scenarios
- Validates theory against exact recomputation
- Mean error: 1.76% for radius perturbations

### Validation

The framework includes **9 comprehensive tests** covering:
- Dispersion curve monotonicity (normal vs anomalous)
- V-parameter scaling relationships
- Mode normalization (∫|Ψ|² dA = 1)
- Physical constraints (n₂ < nₑff < n₁)
- Perturbation accuracy across parameter ranges

**Result:** All 9/9 tests passing ✅

---

## Results & Figures

### 1. Mode Profile

![Mode Profile](figures/02_mode_profile.pdf)
*Figure 2: LP₀₁ electric field amplitude and power density. The field decays exponentially in the cladding (K₀ Bessel function).*

### 2. Dispersion Curve

![Dispersion](figures/03_dispersion_curve.pdf)
*Figure 3: Effective index vs wavelength showing normal material dispersion. The V-parameter decreases with wavelength; mode cutoff occurs at V = 2.405.*

### 3. Perturbation Validation

![Validation](figures/04_radius_validation.pdf)
*Figure 4: Comparison of perturbation theory (blue) vs exact recomputation (red) for core radius variations. Mean error: 1.76%, validating first-order approximation.*

### 4. Index Sensitivity Map

![Sensitivity](figures/05_index_sensitivity.pdf)
*Figure 5: 2D parameter space showing how effective index shifts depend on perturbation magnitude (Δn) and spatial extent (shell width). Contours show regions of constant Δnₑff.*

### 5. Absorption Analysis

![Absorption](figures/06_absorption.pdf)
*Figure 6: Attenuation vs conductivity spanning 6 orders of magnitude. Right panel shows characteristic attenuation length (1/α), relevant for fiber sensor design.*

---

## Key Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Dispersion solver | <10⁻¹⁴ residual | Excellent numerical convergence |
| Mode normalization | 0.21% error | High-quality field profiles |
| Perturbation accuracy | 1.76% mean error | Valid first-order approximation |
| Test coverage | 9/9 passing | Comprehensive validation |
| Code base | 1,922 lines | Production-ready implementation |

---

## Applications

This framework enables:

1. **Fiber sensor design** — Quantify sensitivity to environmental perturbations
2. **Fabrication tolerance analysis** — Predict performance degradation from manufacturing errors
3. **Waveguide optimization** — Guide design choices for integrated photonic circuits
4. **Material characterization** — Extract absorption coefficients from measured attenuation

---

## Skills & Technologies

**Mathematical foundations:**
- Electromagnetic wave theory
- Sturm–Liouville eigenvalue problems
- Bessel functions
- Perturbation theory

**Numerical methods:**
- Root finding (Brent's method)
- Numerical integration
- Eigenvalue solvers

**Software:**
- Python (NumPy, SciPy, Matplotlib)
- LaTeX for technical documentation
- Git version control

---

## What I Learned

This project deepened my understanding of:
- How theoretical physics connects to practical engineering problems
- The power of perturbation methods for sensitivity analysis
- Numerical accuracy and validation strategies
- Creating reproducible research with comprehensive documentation

Most importantly, it demonstrated my ability to **independently tackle complex problems** from theory through implementation to validation — exactly the skills needed for graduate research in photonics and computational electromagnetics.

---

## Explore the Full Project

- **GitHub Repository:** [mmaharebi/fiber-perturbation](https://github.com/mmaharebi/fiber-perturbation)
- **Full Technical Report (PDF):** 37 pages with complete derivations and figures
- **Run the Code:** All results reproducible with `bash run_all.sh`

Feel free to explore the code, run the tests, or reach out if you have questions about the methodology or implementation!

---

## About Me

I'm Mahdy Mahareb, currently preparing applications for funded master's programs in photonics and computational electromagnetics at the University of Kassel, Germany. This project represents my independent work in numerical electromagnetics and demonstrates readiness for graduate-level research.

**Contact:** contact@mahdymahareb.de  
**Website:** https://mahdymahareb.de  
**GitHub:** https://github.com/mmaharebi

---

*Generated: December 2025*  
*Technologies: Python, SciPy, NumPy, Matplotlib, LaTeX*
