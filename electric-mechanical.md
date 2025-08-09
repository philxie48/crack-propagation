# FEM Electro-Mechanical Crack Propagation Simulation Workflow

## Overview
This document outlines the complete workflow for a fully coupled electro-mechanical-phase field crack propagation simulation. The simulation uses a 2-step solution algorithm at each time step:
1. **Electric potential** (sequential)
2. **Mechanical displacement and phase field** (coupled system)

This approach ensures consistent energy balance between mechanical and electrical effects driving crack propagation.

## 1. Material Properties

### Mechanical Properties
- Young's modulus: E = 210 GPa
- Poisson's ratio: ν = 0.3
- Density: ρ = 7380 kg/m³
- **Lamé parameters**: λ = 1.2115×10¹¹ Pa, μ = 8.077×10¹⁰ Pa
- **Body force**: b = [0, -72.4×10³] N/m³ (downward)

### Phase Field Properties
- Critical energy release rate: Gc = 4000 N/m
- Length scale parameter: l₀ = 12.5 μm
- Residual stiffness: k = 1×10⁻⁶

### Electrical Properties
- Electrical conductivity: σ₀ = 8.7×10⁶ S/m
- Vacuum permittivity: ε₀ = 8.854×10⁻¹² F/m
- Constant B: B = 1.24×10⁷ (dimensionless)
- Critical energy release rate (electrical): Gcr^s = 4000 N/m

### Penalty Method Parameters
- Penalty parameter: η = 4×10⁶
- Penalty exponent: m = 2

## 2. Geometry and Mesh
- Domain: 1.0 × 1.0 mm square
- Initial crack: 0.5 mm length × 8 μm width at center (0.5, 0.5) mm
- Crack represented as void in mesh geometry
- High-quality mesh around crack tip: 2.5 μm element size in 0.1 mm radius
- Bulk mesh size: 12.5 μm

## 3. Governing Equations

### 3.1 Electric Potential Equation
The electric potential V is governed by:
```
∇ · (σ(φ) ∇V) = 0
```

Where the degraded electrical conductivity is:
```
σ(φ) = h(φ) σ₀
```

With electrical degradation function:
```
h(φ) = (1 - φ)² + h̄
```
where h̄ is a small residual constant (h̄ = 1×10⁻⁶).

### 3.2 Electric Energy Density
The electric energy density is:
```
W(E) = ½ ε₀ E²
```
where E = -∇V is the electric field strength.

### 3.3 Modified Phase Field Residual
The phase field residual is modified to include electrical effects:

```
Rₑᶠ = ∫_Ω {Gc l₀ (Bₜᶠ)ᵀ ∇φ + [Gc/l₀ + 2ψ(ε) + 2W(E)] Nφ - 2[ψ(ε) + W(E)] - η/(m∆t) ⟨φ̇⟩₋ᵐ} dΩ = 0
```

**Penalty function definition**:
```
⟨x⟩₋ = {-x, x < 0
        {0,  x ≥ 0
```

**Penalty term**:
```
P(φ̇) = η/(m∆t) ⟨φ̇⟩₋ᵐ
```

### 3.4 Critical Energy Release Rate
The critical energy release rate is modified by electrical effects:
```
Gcr^V = Gcr^s - BV²
```
where V is the electric potential.

### 3.5 Mechanical Equilibrium
The mechanical equilibrium remains:
```
∇ · σ = 0
```
with degraded stress:
```
σ = g(φ) D : ε
```

## 4. Finite Element Weak Forms

### 4.1 Electric Potential Weak Form
Find V ∈ H¹ such that:
```
∫_Ω σ(φ) ∇V · ∇δV dΩ = ∫_Γₙᵉ j_n δV dΓ
```
where j_n is the prescribed current density on Neumann boundaries.

### 4.2 Mechanical Equilibrium Residual Form
Find u ∈ [H¹(Ω)]² such that for all test functions δu ∈ [H¹₀(Ω)]²:

```
Rᵉᵘ = ∫_Ω [(1-φ)² + k](Bᵢᵘ)ᵀ σᵢⱼ dΩ - ∫_Ω Nᵀ bⱼ dΩ = 0
```

Where:
- g(φ) = (1-φ)² + k = (1-φ)² + 10⁻⁶
- σᵢⱼ = λtr(ε)δᵢⱼ + 2μεᵢⱼ
- εᵢⱼ = ½(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ)
- bⱼ = [0, -ρg] = [0, -76,518] N/m³ (downward body force)

### 4.3 Phase Field Weak Form
Find φ ∈ H¹ such that:
```
∫_Ω [Gc l₀ ∇δφ · ∇φ + (Gc/l₀ + 2[ψ(ε) + W(E)]) δφ φ + η/(m∆t) δφ ⟨φ̇⟩₋ᵐ] dΩ = 
∫_Ω 2[ψ(ε) + W(E)] δφ dΩ
```

## 5. Boundary Conditions

### 5.1 Mechanical Boundary Conditions
- **Bottom edge (y = 0)**: u = 0 (fixed displacement)
- **Top edge (y = 1 mm)**: ux = 1×10⁻⁵ × t mm, uy = 0 (prescribed shear)
- **Left edge (x = 0)**: Stress-free (natural BC) + one rigid body constraint
- **Right edge (x = 1 mm)**: Stress-free (natural BC) - NO ARTIFICIAL CONSTRAINTS

### 5.1.1 Rigid Body Constraints
- Fix one left boundary node's x-displacement to prevent rigid body translation
- **Important**: No additional constraints on right boundary to prevent artificial stress concentration

### 5.2 Electrical Boundary Conditions
Based on the geometry shown:
- **Left edge (above crack)**: Current injection j_n = j₀ (red region)
- **Left edge (below crack)**: Current extraction j_n = -j₀ (blue region)
- **Right edge**: V = 0 (ground)
- **Top/Bottom edges**: ∂V/∂n = 0 (insulated)
- **Crack surfaces**: ∂V/∂n = 0 (insulated)

### 5.3 Phase Field Boundary Conditions
- **All external boundaries**: ∂φ/∂n = 0 (natural boundary condition)
- **Crack surfaces**: ∂φ/∂n = 0 (natural boundary condition)

## 6. Initial Conditions
- **Electric potential**: V = 0 everywhere
- **Displacement**: u = 0 everywhere  
- **Phase field**: φ = 0 everywhere (crack void represented in mesh geometry)

## 7. Solution Algorithm

### 7.1 Time Stepping
For each time step n+1:

#### Step 1: Solve Electric Potential
1. Assemble electric conductivity matrix: K_V = ∫_Ω σ(φⁿ) ∇N · ∇N dΩ
2. Apply electrical boundary conditions
3. Solve: K_V V^(n+1) = F_V
4. Compute electric field: E^(n+1) = -∇V^(n+1)
5. Compute electric energy density: W(E^(n+1))

#### Step 2: Solve Coupled Mechanical-Phase Field Problem
**Implementation**: The mechanical displacement and phase field are solved together as a coupled system, driven by both mechanical and electrical energy contributions.

1. **Assemble coupled system matrix**:
   ```
   [K_uu  K_uφ] [Δu] = [R_u]
   [K_φu  K_φφ] [Δφ]   [R_φ]
   ```
   Where:
   - K_uu = ∫_Ω g(φⁿ) Bᵀ D B dΩ (mechanical stiffness)
   - K_φφ = ∫_Ω [Gc l₀ ∇N·∇N + (Gc/l₀ + 2[ψ + W]) N·N + η/(m∆t) N·N] dΩ (phase field)
   - K_uφ, K_φu = coupling terms between mechanics and phase field

2. **Apply boundary conditions**:
   - Mechanical: Bottom fixed, top prescribed shear, left/right stress-free
   - Phase field: Natural boundary conditions ∇φ·n = 0

3. **Solve coupled system**: [K_coupled] {u^(n+1), φ^(n+1)} = {F_coupled}

4. **Post-process**:
   - Compute strain: ε^(n+1) = ∇ˢu^(n+1)
   - Compute elastic energy density: ψ(ε^(n+1))
   - Enforce irreversibility: φ^(n+1) = max(φⁿ, φ^(n+1))
   - Enforce bounds: φ^(n+1) ∈ [0,1]

**Note**: This coupled approach ensures consistent energy balance between mechanical and phase field responses to electrical loading.

### 7.2 Solution Strategy
- **Electric potential**: Direct solve of linear conductivity system
- **Mechanical-Phase field**: Direct stiffness matrix approach (no iteration required)
- **Linearization**: Uses current state values for nonlinear terms
- **Robustness**: Single matrix solve per time step, always converges

## 8. Matrix Assembly Details

### 8.1 Electric Potential Matrix
```
K_V[i,j] = ∫_Ω σ(φ) ∇N_i · ∇N_j dΩ
F_V[i] = ∫_Γₙᵉ j_n N_i dΓ
```

### 8.2 Mechanical Matrix
```
K_u[2i,2j] = K_u[2i+1,2j+1] = ∫_Ω g(φ) Bᵢᵀ D B_j dΩ
F_u[2i] = ∫_Γₙᵘ t_x N_i dΓ + ∫_Ω f_x N_i dΩ
F_u[2i+1] = ∫_Γₙᵘ t_y N_i dΓ + ∫_Ω f_y N_i dΩ
```

### 8.3 Phase Field Matrix
```
K_φ[i,j] = ∫_Ω [Gc l₀ ∇N_i · ∇N_j + (Gc/l₀ + 2[ψ + W] + η/(m∆t)) N_i N_j] dΩ
F_φ[i] = ∫_Ω 2[ψ + W] N_i dΩ + η/(m∆t) ∫_Ω (φ̇)^(m-1) N_i dΩ
```

### 8.4 Constitutive Matrix (Plane Strain)
For plane strain linear elastic material:
```
D = E(1-ν)/[(1+ν)(1-2ν)] [  1      ν/(1-ν)     0    ]
                          [ν/(1-ν)     1        0    ]
                          [  0         0    (1-2ν)/2(1-ν)]
```
Numerical values (E = 210 GPa, ν = 0.3):
```
D = 2.538×10¹¹ [1.0  0.429   0  ]
               [0.429  1.0   0  ] Pa
               [ 0     0    0.286]
```

### 8.5 B-Matrix Definitions
**Strain-displacement matrix**:
```
Bᵢᵘ = [∂Nᵢ/∂x    0   ]
      [   0    ∂Nᵢ/∂y]  (3×6 for 3 nodes, 2 DOF each)
      [∂Nᵢ/∂y  ∂Nᵢ/∂x]
```

**Phase field gradient matrix**:
```
Bᵢᶠ = [∂Nᵢ/∂x]
      [∂Nᵢ/∂y]  (2×3 for 3 nodes)
```

**Electric field gradient matrix**:
```
Bᵢᵛ = [∂Nᵢ/∂x]
      [∂Nᵢ/∂y]  (2×3 for 3 nodes)
```

## 9. Implementation Notes

### 9.1 Element Integration
- Use 2×2 Gauss quadrature for triangular elements
- Evaluate all field quantities (σ, ψ, W, g, h) at integration points
- Use nodal interpolation for field variables (V, u, φ)

### 9.2 Degradation Functions
```
g(φ) = (1-φ)² + k    (mechanical degradation)
h(φ) = (1-φ)² + h̄    (electrical degradation)
```

### 9.3 Irreversibility Constraint
Implement using penalty method:
```
R_penalty = η/(m∆t) ∫_Ω N (φ̇)^m dΩ  for φ̇ < 0
```

## 10. Output and Visualization

### 10.1 Field Variables to Monitor
- Electric potential V(x,y,t)
- Electric field magnitude |E(x,y,t)|
- Current density j(x,y,t) = σ(φ)E
- Displacement magnitude |u(x,y,t)|
- Phase field φ(x,y,t)
- Mechanical strain energy density ψ(x,y,t)
- Electric energy density W(x,y,t)
- Total energy density ψ + W

### 10.2 Global Quantities
- Total electric power: P = ∫_Γₙᵉ j_n V dΓ
- Total mechanical work: W_mech = ∫_Γₙᵘ t · u̇ dΓ
- Crack area: A_crack = ∫_Ω φ dΩ
- Total energy dissipation rate

### 10.3 Crack Propagation Metrics
- Crack tip position
- Crack propagation velocity
- Energy release rate at crack tip
- Electrical vs mechanical driving forces

## 11. Computational Parameters
- **Time step**: Δt = 1.0 s (quasi-static loading)
- **Total simulation time**: T = 2000 s (adjustable for electro-mechanical)
- **Convergence tolerance**: 1×10⁻⁶ (relative)
- **Maximum iterations per time step**: 50
- **Integration**: 3-point Gauss quadrature for triangular elements
- **Current density magnitude**: j₀ = 1×10⁸ A/m² (= 100 A/mm²)

## 12. Implementation Requirements
### 12.1 Mesh Specifications
- **Element type**: 3-node linear triangles
- **Crack tip refinement**: 2.5 μm elements within 0.1 mm radius
- **Global mesh size**: 12.5 μm typical element size
- **Quality requirements**: Minimum angle > 15°, aspect ratio < 5

### 12.2 Algorithm Details
- **2-step solution**: Electric → (Mechanical + Phase field coupled) per time step  
- **Electric solver**: Direct sparse solver for conductivity matrix
- **Coupled solver**: Direct stiffness matrix approach (linearized, single solve)
- **No iterations**: Robust linearization eliminates need for convergence loops
- **Matrix format**: Sparse CSR for computational efficiency
- **Coupling treatment**: 
  - Electric → Mechanical: Explicit (electrical energy W drives phase field)
  - Mechanical ↔ Phase field: Implicit (fully coupled within step 2)

This workflow provides a complete framework for implementing the fully coupled electro-mechanical-phase field crack propagation simulation with the electrical effects properly incorporated into the phase field evolution equation.