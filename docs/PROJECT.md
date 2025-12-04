# Perturbation Theory for Step-Index Optical Fibers  
### Eigenvalue Sensitivity to Fabrication Errors via Sturm–Liouville Analysis

## 1. Overview

This project presents a mathematically rigorous, simulation-assisted study of **guided modes in a step-index optical fiber** and their sensitivity to small fabrication-related perturbations.  
The analysis combines:

- **Analytical eigenmode theory** for step-index fibers  
- **Bessel-function dispersion relations**  
- **Sturm–Liouville operator theory**  
- **First-order perturbation methods**  
- **Numerical root-finding and integration**  

The goal is to compute how small deviations in material or geometry—such as **refractive index diffusion**, **core radius errors**, or **weak absorption**—affect the **propagation constant** $\beta$ and **effective index** $n_{\text{eff}} = \beta/k_0$ of the fundamental LP₀₁ mode.

This project is designed as a rigorous and compact demonstration of electromagnetic theory, analytical modeling, and numerical computation relevant to **integrated photonics**, **optical fibers**, **computational electromagnetics**, and **microwave waveguides**.

---

## 2. Physical Problem

We consider a **cylindrical step-index fiber** with:

- Core refractive index: $n_1$  
- Cladding refractive index: $n_2$  
- Core radius: $a$  
- Free-space wavelength: $\lambda_0$

The refractive index profile is:

$$
n(r) =
\begin{cases}
n_1, & r < a,  \\
n_2, & r > a.
\end{cases}
$$

Under the weakly guiding approximation ($n_1 \approx n_2$), the electric field envelope of the LP modes satisfies the scalar Helmholtz equation:

$$
\left[
\frac{1}{r}\frac{d}{dr}\left( r\frac{d}{dr} \right)
+ \left( k_0^2 n^2(r) - \beta^2 \right)
- \frac{m^2}{r^2}
\right] \Psi(r) = 0.
$$

For the fundamental LP₀₁ mode ($m = 0$):

- **Core region** $(r < a)$: $\Psi(r) = A J_0(ur/a)$  
- **Cladding region** $(r > a)$: $\Psi(r) = C K_0(wr/a)$

where

$$
u^2 = (k_0 n_1)^2 - \beta^2, \qquad
w^2 = \beta^2 - (k_0 n_2)^2, \qquad
u^2 + w^2 = V^2.
$$

The **dispersion relation** arises from continuity of $\Psi$ and $\Psi'$ at $r=a$:

$$
\frac{J_0(u)}{u J_1(u)} 
= - \frac{K_0(w)}{w K_1(w)}.
$$

Solving this transcendental equation yields $\beta$ and the effective index $n_{\text{eff}} = \beta/k_0$.

---

## 3. Sturm–Liouville Formulation

The radial wave equation may be written as a **Sturm–Liouville eigenvalue problem**:

$$
\mathcal{L}\Psi
=
-\frac{1}{r}\frac{d}{dr}\left(r\frac{d\Psi}{dr}\right)
+ \left[ \frac{m^2}{r^2} + k_0^2 n^2(r) \right] \Psi
= \beta^2 \Psi.
$$

On the domain $r \in (0, \infty)$, guided modes satisfy:

- Square integrability  
- Exponential decay as $r \to \infty$  
- Vanishing boundary term in Green’s identity  

Thus $\mathcal{L}$ is **self-adjoint**, guaranteeing:

- Real, discrete eigenvalues $\beta^2$ for bound modes  
- Orthogonality of eigenfunctions  
- Validity of **first-order eigenvalue perturbation theory**

This theoretical foundation is rarely shown explicitly in photonics literature, and is a core strength of this project.

---

## 4. First-Order Perturbation Theory

Consider a small perturbation in material or geometry described by:

$$
n^2(r) \rightarrow n^2(r) + \Delta n^2(r), \qquad |\Delta n^2| \ll n^2.
$$

First-order perturbation theory gives the shift in propagation constant:

$$
\Delta \beta
= \frac{\omega}{4 P_z}
\int_0^\infty 
\Delta\epsilon(r)\,|\Psi(r)|^2 \, 2\pi r\, dr,
$$

where $\Delta\epsilon = \epsilon_0\,\Delta n^2$ and $P_z$ is modal power.  
Under scalar power normalization—used widely for LP modes—the denominator simplifies to 1.

Using $\beta = k_0 n_{\text{eff}}$, the shift in effective index is:

$$
\Delta n_{\text{eff}} = \frac{\Delta \beta}{k_0}.
$$

This formula allows sensitivity evaluation for:

- Index diffusion  
- Radius variation  
- Weak absorption (loss)  

---

## 5. Perturbation Scenarios Modeled

### **5.1. Localized Index Perturbation (Doping / Diffusion)**

A diffusion-induced variation near the core boundary:

$$
\Delta n^2(r) =
\begin{cases}
\Delta n_0^2 & a-\delta < r < a+\delta, \\
0 & \text{otherwise},
\end{cases}
$$

is used to compute:

$$
\Delta n_{\text{eff}}(\Delta n_0, \delta).
$$

This models realistic fabrication imperfections in glass or silicon photonics.

---

### **5.2. Core Radius Error ($a \to a + \Delta a$)**

A geometric deviation is approximated via:

1. Numerical derivative:
  $$
  \frac{d\beta}{da} \approx \frac{\beta(a+\Delta a) - \beta(a)}{\Delta a}
  $$
2. Equivalent perturbation ring:
  $$
  \Delta n^2(r) \propto \delta(r-a)\,\Delta a.
  $$

This allows full sensitivity mapping:

$$
n_{\text{eff}}(a+\Delta a) - n_{\text{eff}}(a).
$$

---

### **5.3. Weak Absorption (Complex Index in Core)**

A small conductivity $\sigma$ in the core modifies permittivity:

$$
\Delta\epsilon = -j \frac{\sigma}{\omega}.
$$

Perturbation theory yields the attenuation constant:

$$
\alpha = -2\,\operatorname{Im}(\beta) 
= -2\,\operatorname{Im}(\Delta\beta).
$$

This models **material loss**, **doping absorption**, and **chemical sensing**.

---

## 6. Numerical Implementation

The computational workflow:

### **6.1. Solve unperturbed eigenvalue problem**

- Use SciPy special functions `jv`, `kv`  
- Use bracketing + Brent root-finding  
- Compute $\beta$, $n_{\text{eff}}$, and the LP₀₁ mode profile

### **6.2. Normalize the mode**

$$
\int_0^\infty |\Psi(r)|^2\,2\pi r\,dr = 1.
$$

### **6.3. Evaluate perturbation integrals**

General function:

$$
I = \int_0^\infty \Delta n^2(r)\,|\Psi(r)|^2\,(2\pi r)\,dr.
$$

### **6.4. Compute $\Delta\beta$ and $\Delta n_\text{eff}$**

Use perturbation formula:

$$
\Delta\beta = -\frac{k_0^2}{2\beta} I.
$$

### **6.5. Validation**

For each perturbation:

- Compute $\Delta\beta$ from perturbation theory  
- Modify the physical parameter ($n_1$ or $a$)  
- Re-solve the dispersion equation  
- Compare “exact” vs first-order result

This confirms the validity range of perturbation.

---

## 7. Results Produced

The project generates:

- Effective index $n_{\text{eff}}(\lambda_0)$ curves  
- LP₀₁ mode profile plots  
- Sensitivity maps:
  - $\Delta n_{\text{eff}}$ vs $\Delta n$  
  - $\Delta n_{\text{eff}}$ vs $\Delta a$  
  - Attenuation $\alpha$ vs conductivity $\sigma$  
- Relative error plots comparing perturbation theory with exact recomputation  

These plots demonstrate **analytical correctness**, **numerical stability**, and **real-world relevance**.

---

## 8. Applications

This framework is directly relevant to:

- Optical fiber sensing  
- Integrated photonics process variations (etch & diffusion)  
- High-precision interferometry  
- Mode filtering & modal dispersion engineering  
- Microwave dielectric rod waveguides (mathematically analogous)  

The methodology also generalizes to:

- Graded-index fibers  
- Weakly guiding planar waveguides  
- Perturbative analysis of photonic integrated circuits