# Perturbation Theory for Step-Index Optical Fibers  
### Eigenvalue Sensitivity to Fabrication Errors via Sturm–Liouville Analysis

This repository contains a complete analytical–numerical study of guided modes in a **step-index optical fiber** and their sensitivity to small **fabrication-induced perturbations** in refractive index, geometry, and material absorption.

The project demonstrates:

- Analytical eigenmode theory (LP modes, Bessel & modified Bessel functions)
- Sturm–Liouville formulation of the fiber wave equation
- First-order eigenvalue perturbation theory  
- Numerical solution of transcendental dispersion relations
- Sensitivity of propagation constant and effective index
- Validation through full recomputation of eigenvalues

It is designed as a rigorous, compact demonstration of electromagnetic modeling suitable for **integrated photonics**, **fiber optics**, **computational electromagnetics**, and **microwave waveguides**.

---

## ✨ Key Features

- **LP$_{01}$ Mode Solver**  
  Computes $\beta$ and the effective index by solving the classic step-index dispersion relation using SciPy's Bessel functions.

- **Mode Profile Generation**  
  Produces normalized radial field profiles using $J_0$ and $K_0$ solutions.

- **First-Order Perturbation Engine**  
  Implements the Sturm–Liouville eigenvalue perturbation formula:
  $$
  \Delta\beta = (\omega / 4P_z) \int \Delta\epsilon(r)|\Psi(r)|^2 2\pi r dr
  $$
  
- **Three Perturbation Scenarios**
  - Localized index change (doping / diffusion)
  - Core radius error ($a \to a + \Delta a$)
  - Small absorption in the core (computes attenuation $\alpha$)

- **Validation Tools**  
  Perturbation results are compared against recomputation of β with modified parameters.

---

## 📁 Repository Structure

```
fiber-perturbation/
├── PROJECT.md
├── README.md
├── docs/
│   └── theory.md
├── src/
│   ├── fiber_modes.py
│   ├── perturbation.py
│   └── utils.py
├── notebooks/
│   ├── 01_unperturbed_modes.ipynb
│   ├── 02_index_perturbation.ipynb
│   ├── 03_radius_error.ipynb
│   └── 04_loss_perturbation.ipynb
└── figures/
    ├── lp01_profile.png
    ├── neff_vs_lambda.png
    ├── delta_neff_vs_dn.png
    ├── delta_neff_vs_da.png
    ├── attenuation_vs_sigma.png
    └── perturbation_validation.png
```
