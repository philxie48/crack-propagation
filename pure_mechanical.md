# Mechanical Crack Propagation: FEM Weak Form with Penalty Method

## 1. Parameters

### Material Properties
- **Young's modulus**: E = 210 GPa, **Poisson's ratio**: ν = 0.3, **Density**: ρ = 7380 kg/m³
- **Lamé parameters**: λ = 1.2115×10¹¹ Pa, μ = 8.077×10¹⁰ Pa
- **Body force**: b = [0, -72.4×10³] N/m³ (downward)

### Phase Field Properties
- **Length scale**: l₀ = 1.25×10⁻⁵ m (12.5 μm)
- **Critical energy release rate**: Gc = 4000 N/m
- **Degradation function**: g(φ) = (1-φ)² + k, k = 10⁻⁶

### Penalty Method Parameters (for Irreversibility)
- **Penalty parameter**: η = 4×10⁶
- **Penalty exponent**: m = 2
- **Time step**: Δt = 1.0 s (quasi-static loading)

### Numerical Parameters
- **Convergence tolerance**: 1×10⁻⁶ (relative)
- **Maximum iterations**: 50 per time step
- **Integration**: 3-point Gauss quadrature for triangular elements
- **Total simulation time**: T = 100 s (adjustable)

### Geometry
- **Domain**: 1.0×1.0 mm square
- **Initial crack**: 0.5 mm length × 8 μm width at center (0.5, 0.5) mm
- **Crack representation**: Void in mesh geometry (φ = 0 in crack region)

## 2. Initialization

- **Phase field**: φ = 0 everywhere (intact material) 
  - **Note**: Crack void is represented directly in mesh geometry, not through φ = 1
  - Phase field φ tracks crack propagation beyond the initial void
- **Displacement**: u = 0 everywhere

## 3. Boundary Conditions

### Mechanical Field
- **Bottom edge (y = 0)**: u = 0 (fixed displacement)
- **Top edge (y = 1 mm)**: uₓ = 1×10⁻⁵t mm, uᵧ = 0 (prescribed shear)
- **Left edge (x = 0)**: stress-free (natural BC) + one rigid body constraint
- **Right edge (x = 1 mm)**: stress-free (natural BC) - NO ARTIFICIAL CONSTRAINTS

### Rigid Body Constraints
- Fix one left boundary node's x-displacement to prevent rigid body translation
- **Important**: No additional constraints on right boundary to prevent artificial stress concentration

### Phase Field
- **All boundaries**: ∇φ·n = 0 (natural BC)

## 4. FEM Weak Form Formulation

### 4.1 Mechanical Equilibrium Residual Form

Find u ∈ [H¹(Ω)]² such that for all test functions δu ∈ [H¹₀(Ω)]²:

```
Rᵤ(u, δu) = ∫_Ω g(φ)σᵢⱼ(u):∇δu dΩ - ∫_Ω bⱼ·δu dΩ = 0
```

Where:
- g(φ) = (1-φ)² + k = (1-φ)² + 10⁻⁶
- σᵢⱼ = λtr(ε)δᵢⱼ + 2μεᵢⱼ
- εᵢⱼ = ½(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ)
- bⱼ = [0, -72.4×10³] N/m³

**Substituted form**:
```
Rᵤ = ∫_Ω [(1-φ)² + 10⁻⁶][1.2115×10¹¹tr(ε)I + 2×8.077×10¹⁰ε]:∇δu dΩ 
   - ∫_Ω [0, -72.4×10³]·δu dΩ = 0
```

### 4.2 Phase Field Residual Form with Penalty Method

Based on the penalty method approach, the coupled system becomes:

Find {u, φ} such that for all test functions {δu, δφ}:

```
Rᵉᵘ = ∫_Ω [(1-φ)² + k](Bᵢᵘ)ᵀ σᵢⱼ dΩ - ∫_Ω Nᵀ bⱼ dΩ =0
```

```
Rᵉᶠ = ∫_Ω {Gc l₀(Bᵢᶠ)ᵀ ∇φ + [Gc/l₀ + 2ψ(ε)] Nφ - 2ψ(ε) - η/(mΔt) ⟨φ̇⟩₋ᵐ} dΩ =0
```

**Penalty function definition**:
```
⟨x⟩₋ = {-x, x < 0
         {0,  x ≥ 0
```

**Penalty term**:
```
P(φ̇) = η/(mΔt) ⟨φ̇⟩₋ᵐ
```

Where:
- **η** = 4×10⁶ (penalty parameter)
- **m** = 2 (penalty exponent)
- **Δt** = time step
- **ψ(ε)** = [μ(ε:ε) + (1/2)λ(tr ε)²]⁺ (positive elastic energy density)
- **Gc** = 4000 N/m (fracture energy)
- **l₀** = 1.25×10⁻⁵ m (length scale)

### B-Matrix Definitions

For 2D plane strain problems:

**Phase field gradient matrix**:
```
Bᵢᶠ = [∂Nᵢ/∂x]
      [∂Nᵢ/∂y]
```

**Strain-displacement matrix**:
```
Bᵢᵘ = [∂Nᵢ/∂x    0   ]
      [   0    ∂Nᵢ/∂y]
      [∂Nᵢ/∂y  ∂Nᵢ/∂x]
```

### Constitutive Matrix D

For plane strain linear elastic material:
```
D = E(1-ν)/[(1+ν)(1+2ν)] [  1      ν/(1-ν)     0    ]
                          [ν/(1-ν)     1        0    ]
                          [  0         0    (1-2ν)/2(1-ν)]
```

Substituting material properties (E = 210 GPa, ν = 0.3):
```
D = 2.538×10¹¹ [1.0  0.429   0  ]
               [0.429  1.0   0  ]
               [ 0     0    0.286]
```

## 5. FEM Solution Workflow

### Simultaneous Coupled Solution Strategy

The system is solved as a **fully coupled linear problem** using direct stiffness matrix approach. Both displacement and phase field are solved simultaneously at each time step using linearization.

### Direct Stiffness Matrix Solution Algorithm
1. **Initialize**: u = 0, φ = 0 everywhere (Section 2)
2. **For each time step t**:
   
   a. **Assemble coupled system matrix**:
      - Compute stiffness matrix K and force vector F
      - Mechanical block: K_uu = ∫ g(φ) B^T D B dΩ
      - Phase field block: K_φφ = ∫ [Gc l₀ ∇N·∇N + (Gc/l₀ + 2ψ) N·N] dΩ - Penalty
      - Coupling blocks: K_uφ and K_φu from degradation and energy derivatives
      - Force vector: F_u (body forces), F_φ (energy driving forces + penalty)
   
   b. **Apply boundary conditions**:
      - Dirichlet: Set rows/columns, modify F
      - Rigid body constraints to prevent spurious motion
   
   c. **Solve directly**: K x = F using sparse solver
   
   d. **Update solution**:
      - Extract u and φ from solution vector x
      - Apply bounds: φ ∈ [0,1]
      - Enforce irreversibility: φ_new = max(φ_old, φ_new)
   
   e. **Time stepping**: Update φ_old ← φ, advance t ← t + Δt

### Key Implementation Notes
- **Linearized approach**: Treat as quasi-linear problem at each time step
- **Direct solution**: Single matrix solve per time step (no iteration)
- **Proven robustness**: Based on successful electro-mechanical implementation
- **Initialization**: Simple u = 0, φ = 0 is sufficient
- **Boundary conditions**: Applied in standard matrix form K x = F

## 6. Matrix Assembly Details

### Coupled System Matrix (from Penalty Method)

The coupled system matrix has the form:
```
[Kᵘᵘ  Kᵘᶠ] [δu] = [Rᵘ]
[Kᶠᵘ  Kᶠᶠ] [δφ]   [Rᶠ]
```

Where the submatrices are derived from:
```
Kᵢⱼᵘᵘ = ∂Rᵢᵘ/∂u,  Kᵢⱼᵘᶠ = ∂Rᵢᵘ/∂φ
Kᵢⱼᶠᵘ = ∂Rᵢᶠ/∂u,  Kᵢⱼᶠᶠ = ∂Rᵢᶠ/∂φ
```

### Mechanical Matrix Components
```
Kᵤᵢⱼᵘᵘ = ∫_Ω g(φ) BᵢᵀD Bⱼ dΩ
Kᵤᵢⱼᵘᶠ = ∫_Ω [-2(1-φ)σᵢⱼ + 2ψ(ε)] BᵢᵀNⱼ dΩ
Fᵤᵢ = ∫_Ω b·Nᵢ dΩ
```

### Phase Field Matrix Components
```
Kφᵢⱼᶠᵘ = ∫_Ω 2∂ψ(ε)/∂ε Nᵢ Bⱼ dΩ
Kφᵢⱼᶠᶠ = ∫_Ω [Gc l₀∇Nᵢ·∇Nⱼ + (Gc/l₀)NᵢNⱼ + 2∂ψ(ε)/∂φ NᵢNⱼ] dΩ - Penalty terms
Fφᵢ = ∫_Ω 2ψ(ε)Nᵢ dΩ
```

### Penalty Term Contribution
```
Ppenaltyᵢⱼ = ∫_Ω η/(mΔt) m⟨φ̇⟩₋^(m-1) ∂⟨φ̇⟩₋/∂φ NᵢNⱼ dΩ
```

Where:
- **Nᵢ**: Shape functions
- **Bᵢ**: Strain-displacement matrix  
- **D**: Constitutive matrix
- **g(φ)** = (1-φ)² + k (degradation function)
- **η** = 4×10⁶, **m** = 2 (penalty parameters)

## 7. Implementation Requirements

### 7.1 Mesh Specifications
- **Element type**: 3-node linear triangles
- **Crack tip refinement**: 2.5 μm elements within 0.1 mm radius
- **Global mesh size**: 12.5 μm typical element size
- **Quality requirements**: Minimum angle > 15°, aspect ratio < 5

### 7.2 Algorithm Details

**Two Implementation Approaches Available:**

#### Approach 1: Direct Stiffness Matrix (Recommended)
- **File**: `fem_crack_simulation_updated.py`
- **Method**: Linearized direct solve (single matrix solve per time step)
- **Robustness**: Always converges, no iteration required
- **Performance**: ~15 seconds per time step for typical problems
- **Basis**: Proven successful from electro-mechanical implementation

#### Approach 2: Newton-Raphson Iteration (Exact Nonlinear)
- **File**: `fem_crack_simulation_newton_raphson.py`
- **Method**: Iterative residual-based nonlinear solver
- **Convergence**: Requires tolerance checking (1×10⁻⁶)
- **Performance**: Variable (depends on convergence rate)
- **Use case**: When exact nonlinear solution is required

**Common Features:**
- **Time integration**: Backward Euler for phase field evolution
- **Matrix format**: Sparse CSR for computational efficiency
- **Boundary condition enforcement**: Direct elimination method

### 7.3 Validation Metrics
- **Energy conservation**: (excluding dissipation)
- **Irreversibility check**: φ̇ ≥ 0 everywhere
- **Physical behavior**: Stress concentration at crack tip
- **Convergence**: Mesh and time step independence

### 7.4 Post-Processing Quantities
- **Strain energy density**: ψ = ½ε:D:ε
- **Stress field**: σ = g(φ)D:ε  
- **Crack area**: A_crack = ∫_Ω φ dΩ
- **Energy release rate**: G = -∂U/∂a (numerical differentiation)
