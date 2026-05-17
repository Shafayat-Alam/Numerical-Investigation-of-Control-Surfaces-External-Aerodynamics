# Numerical Investigation of Rocket Fin Aerodynamics
## CFD Analysis at Mach 0.7 using OpenFOAM v2412

**Shafayat Alam**  
*Undergraduate Researcher*  
Advanced Fluid Dynamics, Propulsion, and Energy Lab (AFPEL)  
Department of Mechanical Engineering, Stony Brook University

---

## Executive Summary

This study compares the aerodynamic performance of two rocket fin geometries at Mach 0.7 using high-fidelity CFD simulation. A blunt-edged baseline fin and a streamlined aerodynamic fin were analyzed using steady-state RANS equations with Spalart-Allmaras turbulence modeling in OpenFOAM v2412.

### Key Results

| Metric | Blunt Fin | Aero Fin | Change |
|--------|-----------|----------|--------|
| **Total Drag ($C_D$)** | 0.02970 | 0.02057 | **-30.7%** |
| **Total Lift ($C_L$)** | 0.19837 | 0.07459 | -62.4% |
| **L/D Ratio** | 6.68 | 3.63 | -45.7% |
| **Pressure Drag** | 0.01769 (59.6%) | 0.01252 (60.9%) | -29.2% |
| **Viscous Drag** | 0.01201 (40.4%) | 0.00804 (39.1%) | -33.1% |

**Physical Insight:** The 30.7% drag reduction is entirely due to elimination of base separation wake. The streamlined fin sacrifices lift due to strong counter-acting viscous forces (-0.10368) on the curved trailing edge.

---

## Problem Definition

### Geometry

**Trapezoidal Fin Specifications:**

| Parameter | Value |
|-----------|-------|
| Root chord ($c_r$) | 5.50 in (0.1397 m) |
| Tip chord ($c_t$) | 2.144 in (0.0545 m) |
| Semi-span ($b$) | 2.75 in (0.0699 m) |
| Thickness ($t$) | 0.125 in (0.00318 m) |
| Aspect ratio (AR) | 0.720 |
| Planform area ($S$) | 10.51 in² (0.006780 m²) |
| Mean aerodynamic chord ($c_{\text{MAC}}$) | 4.067 in (0.1033 m) |

$$c_{\text{MAC}} = \frac{2}{3}\left(c_r + c_t - \frac{c_r c_t}{c_r + c_t}\right)$$

![Fin Geometry](images/fin_planform.png)

### Two Configurations

**Blunt Fin:** Rectangular cross-section with flat 90° leading and trailing edges

**Aero Fin:** Curved leading edge with cusped trailing edge for streamlined pressure recovery

![Configuration Comparison](images/fin_configurations.png)

---

## Flow Conditions

| Parameter | Value |
|-----------|-------|
| Mach number ($M$) | 0.7 |
| Altitude | 2,000 ft (ISA) |
| Freestream velocity ($U_\infty$) | 236.6 m/s |
| Static pressure ($p_\infty$) | 94,212 Pa |
| Static temperature ($T_\infty$) | 283.19 K |
| Density ($\rho_\infty$) | 1.154 kg/m³ |
| Dynamic pressure ($q_\infty$) | 32,298 Pa |
| Reynolds number ($\text{Re}_{\text{MAC}}$) | $1.58 \times 10^6$ |
| Angle of attack | 0° |

**Speed of sound:**

$$a = \sqrt{\gamma RT} = \sqrt{1.4 \times 287.05 \times 283.19} = 337.4 \text{ m/s}$$

**Freestream velocity:**

$$U_\infty = M \cdot a = 0.7 \times 337.4 = 236.2 \text{ m/s}$$

---

## Governing Equations

### Reynolds-Averaged Navier-Stokes (RANS)

All flow variables decomposed into mean and fluctuating components:

$$\phi = \bar{\phi} + \phi'$$

**Continuity (Mass Conservation):**

$$\nabla \cdot (\bar{\rho}\bar{\mathbf{u}}) = 0$$

**Momentum Conservation:**

$$\nabla \cdot (\bar{\rho}\bar{\mathbf{u}} \otimes \bar{\mathbf{u}}) = -\nabla \bar{p} + \nabla \cdot \bar{\boldsymbol{\tau}} - \nabla \cdot (\bar{\rho}\overline{\mathbf{u}' \otimes \mathbf{u}'})$$

Where the Reynolds stress tensor is:

$$\boldsymbol{\tau}^R = -\bar{\rho}\overline{\mathbf{u}' \otimes \mathbf{u}'} = \mu_t[\nabla \bar{\mathbf{u}} + (\nabla \bar{\mathbf{u}})^T] - \frac{2}{3}\bar{\rho}k\mathbf{I}$$

**Energy (Total Enthalpy):**

$$\nabla \cdot (\bar{\rho}\bar{\mathbf{u}}\bar{h}_0) = \nabla \cdot (\alpha_{\text{eff}} \nabla \bar{h}_0)$$

Where $h_0 = h + \frac{1}{2}|\mathbf{u}|^2$ and $\alpha_{\text{eff}} = \alpha + \alpha_t$

### Equation of State

**Perfect Gas Law:**

$$\rho = \frac{p}{RT}$$

$$h = C_p T$$

Where $R = 287.05$ J/kg·K and $C_p = 1005$ J/kg·K

**Temperature recovery:**

$$T = \frac{h_0 - \frac{1}{2}|\mathbf{u}|^2}{C_p}$$

### Viscosity Model

**Sutherland's Law:**

$$\mu(T) = \mu_0 \left(\frac{T}{T_0}\right)^{3/2} \frac{T_0 + S}{T + S}$$

Where $\mu_0 = 1.716 \times 10^{-5}$ kg/m·s, $T_0 = 273.15$ K, $S = 110.4$ K

**Prandtl number:** $\text{Pr} = 0.72$

---

## Turbulence Model: Spalart-Allmaras

Single-equation model solving for modified turbulent viscosity $\tilde{\nu}$:

$$\frac{\partial (\rho \tilde{\nu})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \tilde{\nu}) = \frac{1}{\sigma_{\tilde{\nu}}}\nabla \cdot [\rho(\nu + \tilde{\nu})\nabla \tilde{\nu}] + \frac{c_{b2}}{\sigma_{\tilde{\nu}}}\rho(\nabla \tilde{\nu})^2 + \rho \tilde{P}_{\tilde{\nu}} - \rho \tilde{D}_{\tilde{\nu}}$$

**Production term:**

$$\tilde{P}_{\tilde{\nu}} = c_{b1}(1 - f_{t2})\tilde{S}\tilde{\nu}$$

$$\tilde{S} = S + \frac{\tilde{\nu}}{\kappa^2 d^2}f_{v2}$$

Where $S = \sqrt{2\Omega_{ij}\Omega_{ij}}$ (vorticity magnitude) and $d$ = distance to wall

**Destruction term:**

$$\tilde{D}_{\tilde{\nu}} = \left(c_{w1}f_w - \frac{c_{b1}}{\kappa^2}f_{t2}\right)\left(\frac{\tilde{\nu}}{d}\right)^2$$

**Eddy viscosity:**

$$\nu_t = \tilde{\nu}f_{v1}, \quad f_{v1} = \frac{\chi^3}{\chi^3 + c_{v1}^3}, \quad \chi = \frac{\tilde{\nu}}{\nu}$$

**Model constants:**
- $\sigma_{\tilde{\nu}} = 2/3$
- $c_{b1} = 0.1355$
- $c_{b2} = 0.622$
- $\kappa = 0.41$
- $c_{v1} = 7.1$
- $c_{w1} = c_{b1}/\kappa^2 + (1 + c_{b2})/\sigma_{\tilde{\nu}}$
- $c_{w2} = 0.3$
- $c_{w3} = 2.0$

---

## Numerical Method

### Solver Configuration

**OpenFOAM solver:** `rhoSimpleFoam` (steady-state, compressible)

**Algorithm:** SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)

### Discretization Schemes

**Temporal:**
- Steady-state (no time derivative)

**Gradients:**
- Gauss linear (2nd order central differencing)

**Convection (divergence terms):**
- Velocity, enthalpy: `linearUpwind` (2nd order upwind)
- Turbulence ($\tilde{\nu}$): `upwind` (1st order)
- Density flux: `upwind` (1st order, stability)

**Diffusion (Laplacian terms):**
- `Gauss linear corrected` (2nd order with non-orthogonal correction)

### Linear Solvers

**Pressure ($p$):**
- Solver: GAMG (Geometric-Algebraic Multi-Grid)
- Smoother: GaussSeidel
- Tolerance: $10^{-6}$

**Transport variables ($\mathbf{u}$, $h$, $\tilde{\nu}$):**
- Solver: PBiCGStab (Preconditioned Bi-Conjugate Gradient Stabilized)
- Preconditioner: DILU (Diagonal Incomplete LU)
- Tolerance: $10^{-6}$

### Under-Relaxation Factors

Applied to ensure stability in nonlinear iterations:

- Pressure: 0.2
- Velocity: 0.5
- Density: 0.01
- Enthalpy: 0.5
- Turbulence: 0.5

### Boundary Conditions

| Variable | Inlet | Outlet | Fin Surface | Symmetry |
|----------|-------|--------|-------------|----------|
| Velocity ($\mathbf{u}$) | Fixed (236.6, 0, 0) m/s | Zero gradient | No-slip (0, 0, 0) | Symmetry |
| Pressure ($p$) | Zero gradient | Fixed 94,212 Pa | Zero gradient | Symmetry |
| Temperature ($T$) | Fixed 283.19 K | Zero gradient | Adiabatic wall | Symmetry |
| S-A variable ($\tilde{\nu}$) | Fixed $4.5 \times 10^{-5}$ m²/s | Zero gradient | Wall (0) | Symmetry |
| Eddy viscosity ($\nu_t$) | Calculated | Calculated | Wall function | Symmetry |

**Wall functions used:**
- Turbulent thermal diffusivity: `alphatJayatillekeWallFunction`
- Turbulent viscosity: `nutLowReWallFunction`

---

## Computational Mesh

### Mesh Statistics

| Property | Value |
|----------|-------|
| Total cells | 5,879,477 |
| Total nodes | 6,324,343 |
| Hexahedral cells | 5,595,244 (95.2%) |
| Polyhedral cells | 268,697 (4.6%) |
| Domain volume | 4.798 m³ |

### Mesh Quality Metrics

| Metric | Value | Status |
|--------|-------|--------|
| Max aspect ratio | 58.59 | ✓ (in boundary layer) |
| Max face skewness | 3.98 | ✓ (< 4.0 threshold) |
| Average non-orthogonality | 8.47° | ✓ Excellent |
| Max non-orthogonality | 70.00° | ✓ (at limit) |

### Boundary Layer Resolution

**Dimensionless wall distance:**

$$y^+ = \frac{\rho u_\tau y}{\mu}, \quad u_\tau = \sqrt{\frac{\tau_w}{\rho}}$$

**Aero Fin:**
- Area-weighted average: $y^+ = 2.905$ (viscous sublayer)
- Minimum: $y^+ = 0.049$
- Maximum: $y^+ = 103.95$ (leading edge stagnation)

**Blunt Fin:**
- Area-weighted average: $y^+ = 15.757$ (buffer layer)
- Minimum: $y^+ = 2.217$
- Maximum: $y^+ = 133.86$ (trailing edge corners)

**Physical interpretation:** Aero fin achieves direct viscous sublayer resolution ($y^+ < 5$), enabling accurate skin friction calculation without wall functions.

![Mesh Visualization](images/mesh_visualization.png)

---

## Results

### Aero Fin Performance

**Convergence:**
- Total iterations: 2,498
- Final velocity residual: $6.65 \times 10^{-8}$
- Final pressure residual: $4.78 \times 10^{-5}$

![Aero Fin Residuals](images/aero_residuals.png)

**Force Coefficients:**

$$C_D = \frac{D}{q_\infty S}, \quad C_L = \frac{L}{q_\infty S}$$

| Coefficient | Pressure | Viscous | **Total** |
|-------------|----------|---------|-----------|
| Drag ($C_D$) | 0.01252 | 0.00804 | **0.02057** |
| Lift ($C_L$) | 0.17828 | -0.10368 | **0.07459** |

**Key observation:** Negative viscous lift component indicates adverse pressure gradient on curved surface creates downward-acting shear forces.

![Aero Fin Forces](images/aero_forces.png)

### Blunt Fin Performance

**Convergence:**
- Total iterations: 1,013
- Final velocity residual: $1.63 \times 10^{-8}$
- Final pressure residual: $1.01 \times 10^{-6}$

![Blunt Fin Residuals](images/blunt_residuals.png)

**Force Coefficients:**

| Coefficient | Pressure | Viscous | **Total** |
|-------------|----------|---------|-----------|
| Drag ($C_D$) | 0.01769 | 0.01201 | **0.02970** |
| Lift ($C_L$) | 0.09188 | 0.10648 | **0.19837** |

**Key observation:** Both pressure and viscous components contribute constructively to lift.

![Blunt Fin Forces](images/blunt_forces.png)

---

## Comparative Analysis

### Drag Reduction

**Total drag change:**

$$\Delta C_D = \frac{C_{D,\text{blunt}} - C_{D,\text{aero}}}{C_{D,\text{blunt}}} = \frac{0.02970 - 0.02057}{0.02970} = 30.7\%$$

**Component breakdown:**

| Component | Blunt | Aero | Reduction |
|-----------|-------|------|-----------|
| Pressure drag | 0.01769 | 0.01252 | 29.2% |
| Viscous drag | 0.01201 | 0.00804 | 33.1% |
| **Total drag** | 0.02970 | 0.02057 | **30.7%** |

**Base drag contribution:**

Blunt trailing edge creates low-pressure recirculation zone:

$$C_{D,\text{base}} = \frac{(p_\infty - p_{\text{base}}) A_{\text{base}}}{q_\infty S} \approx 0.0177$$

This accounts for **59.5% of total blunt fin drag**.

### Lift Trade-off

**Total lift change:**

$$\Delta C_L = \frac{C_{L,\text{aero}} - C_{L,\text{blunt}}}{C_{L,\text{blunt}}} = \frac{0.07459 - 0.19837}{0.19837} = -62.4\%$$

**Component breakdown:**

| Component | Blunt | Aero | Change |
|-----------|-------|------|--------|
| Pressure lift | 0.09188 | 0.17828 | +94.1% |
| Viscous lift | 0.10648 | -0.10368 | -197.4% |
| **Total lift** | 0.19837 | 0.07459 | **-62.4%** |

**Physical mechanism:**

For surface element at angle $\theta$:

$$dL = (-p \sin\theta + \tau_w \cos\theta) dA$$

- **Blunt fin:** Small $\theta$ → both terms positive (constructive)
- **Aero fin:** Large $\theta$ on curved surface → $\tau_w$ term dominates and is negative (destructive)

### Efficiency Comparison

**Lift-to-drag ratio:**

$$\frac{L}{D} = \frac{C_L}{C_D}$$

| Configuration | $C_L$ | $C_D$ | $L/D$ |
|---------------|-------|-------|-------|
| Blunt Fin | 0.19837 | 0.02970 | **6.68** |
| Aero Fin | 0.07459 | 0.02057 | **3.63** |

**Design implications:**

- **High-speed cruise (minimize drag):** Aero fin preferred (30.7% drag reduction)
- **Control maneuvers (maximize lift):** Blunt fin preferred (166% more lift)
- **Optimal design:** Mission-dependent; may require variable geometry

---

## Flow Physics Interpretation

### Compressibility Effects

At leading edge stagnation ($u = 0$), compressible isentropic relation:

$$C_{p,\text{stag}} = \frac{2}{\gamma M^2}\left[\left(1 + \frac{\gamma-1}{2}M^2\right)^{\gamma/(\gamma-1)} - 1\right]$$

For $M = 0.7$:

$$C_{p,\text{stag}} = 1.076$$

Incompressible prediction: $C_p = 1.0$

**Result:** 7.6% pressure increase due to compressibility at stagnation point.

### Transonic Pocket Formation

**Critical pressure coefficient** (sonic condition):

$$C_p^* \approx -0.92 \text{ at } M = 0.7$$

If local $C_p < -0.92$, flow becomes locally supersonic.

Sharp leading edge on blunt fin likely creates $C_p \approx -1.5$ → local supersonic pocket forms even at freestream $M = 0.7$.

### Boundary Layer Physics

**Viscous sublayer** ($y^+ < 5$):

$$u^+ = y^+$$

**Log-law layer** ($30 < y^+ < 300$):

$$u^+ = \frac{1}{\kappa}\ln(y^+) + C$$

Where $\kappa = 0.41$, $C \approx 5.0$

**Aero fin:** Resolves viscous sublayer → direct calculation of wall shear stress:

$$\tau_w = \mu\left(\frac{\partial u}{\partial y}\right)_{\text{wall}}$$

---

## Conclusions

### Summary of Findings

1. **Streamlined geometry reduces total drag by 30.7%** through elimination of base separation wake (59.5% of baseline drag)

2. **Streamlined geometry reduces total lift by 62.4%** due to strong counter-acting viscous forces on curved trailing edge

3. **Blunt fin provides higher raw lift** (0.19837 vs 0.07459) with both pressure and viscous components acting constructively

4. **Aero fin provides superior drag efficiency** but at cost of reduced control authority

5. **Design trade-off is mission-dependent:**
   - High-speed flight: Aero fin (minimize drag)
   - Control maneuvers: Blunt fin (maximize lift)

### Design Recommendations

**For transonic rocket applications:**
- Use streamlined fins during high-speed ascent (drag-critical)
- Consider blunt fins for control-intensive phases
- Investigate variable geometry or deployable control surfaces
- Validate with wind tunnel testing at matched Reynolds number

---

## Future Work

### Verification & Validation

**Grid Convergence Study:**

Richardson Extrapolation with three grid levels:

$$\text{GCI} = \frac{F_s |\varepsilon|}{r^p - 1}$$

Where $F_s = 1.25$, $r$ = refinement ratio, $p$ = observed order of accuracy

**Experimental Validation:**
- Wind tunnel testing at $M = 0.7$, $\text{Re} = 1.6 \times 10^6$
- 6-component force balance measurements
- Surface pressure distribution (20-30 pressure taps)
- Comparison criteria: $C_D$, $C_L$ within ±5%

### Shape Optimization

**Global optimization:**
- Evolutionary genetic algorithm for design space exploration
- Objective function: minimize drag OR maximize L/D ratio
- Constraints: minimum lift requirement, structural thickness

**Local refinement:**
- Adjoint-based gradient optimization
- Design variables: leading edge radius, trailing edge cusp angle, thickness distribution
- Computational cost: single adjoint solve = single flow solve

### Sensitivity Analysis

**Global sensitivity (Sobol indices):**

$$S_i = \frac{\text{Var}_i(E_{-i}(C_D | x_i))}{\text{Var}(C_D)}$$

Identifies which geometric parameters dominate aerodynamic performance variance.

**Target parameters:** root chord, taper ratio, thickness, sweep angle, aspect ratio

---

## Computational Resources

**Hardware:**
- HPC Cluster: "Intruder" (Stony Brook University)
- CPU: Intel Xeon (specific model not disclosed)
- RAM: Distributed across compute nodes
- Storage: High-performance parallel filesystem

**Computational Cost:**

| Configuration | Iterations | Wall Time | Core-hours |
|---------------|------------|-----------|------------|
| Blunt Fin | 1,013 | ~8 hours | ~192 |
| Aero Fin | 2,498 | ~18 hours | ~432 |

**Software:**
- OpenFOAM v2412 (open-source)
- ParaView for post-processing visualization
- Python (NumPy, Matplotlib) for data analysis

---

## References

1. **Versteeg, H.K. & Malalasekera, W.** (2007)  
   *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*  
   Pearson Prentice Hall, 2nd Edition

2. **Spalart, P.R. & Allmaras, S.R.** (1994)  
   "A One-Equation Turbulence Model for Aerodynamic Flows"  
   *La Recherche Aérospatiale*, No. 1, pp. 5-21

3. **Moukalled, F., Mangani, L., & Darwish, M.** (2016)  
   *The Finite Volume Method in Computational Fluid Dynamics*  
   Springer International Publishing

4. **NOAA, NASA, & USAF** (1976)  
   "U.S. Standard Atmosphere, 1976"  
   NASA-TM-X-74335

5. **Stern, F., Wilson, R.V., Coleman, H.W., & Paterson, E.G.** (1999)  
   "Verification and Validation of CFD Simulations"  
   IIHR Report No. 407

---

## Acknowledgements

**Primary Advisor:**  
Dr. Foluso Ladeinde (Stony Brook University) — Research guidance and HPC access

**Pedagogical Foundation:**  
Dr. Spencer Zimmerman (Stony Brook University) — Fluid mechanics instruction

**Technical Resources:**
- Dr. Aidan Wimshurst (Fluid Mechanics 101) — OpenFOAM implementation
- H.K. Versteeg & W. Malalasekera — Finite Volume Method fundamentals
- F. Moukalled et al. — Advanced OpenFOAM applications

---

**Author:** Shafayat Alam  
**Institution:** Stony Brook University, Department of Mechanical Engineering  
**Date:** May 2026  
**Version:** 1.0

---
