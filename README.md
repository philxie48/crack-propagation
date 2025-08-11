# Phase Field Crack Propagation Simulation Suite

A comprehensive finite element method (FEM) implementation for simulating crack propagation using phase field modeling, supporting both pure mechanical and electro-mechanical driving forces.

## 📋 Overview

This repository provides two complementary crack propagation simulation approaches:

1. **Pure Mechanical Simulation** - Traditional phase field fracture mechanics
2. **Electro-Mechanical Simulation** - Coupled electrical-mechanical crack propagation

Both simulations use the same robust mesh generation system and implement advanced phase field modeling with penalty-based irreversibility constraints.

## 🏗️ Simulation Workflow

### Step 1: Mesh Generation
```bash
python robust_crack_mesh.py
```

**Purpose**: Generates high-quality finite element meshes with pre-defined crack geometry
- **Output**: `robust_crack.msh` (GMSH format), `robust_crack.npz` (NumPy format)
- **Features**: Adaptive mesh refinement around crack tip, void representation in geometry
- **Crack Geometry**: 0.5 mm length × 8 μm width centered in 1.0×1.0 mm domain

### Step 2A: Pure Mechanical Simulation

**Option 1: Direct Stiffness Matrix (Recommended)**
```bash
python fem_crack_simulation_updated.py
```

**Option 2: Newton-Raphson Iteration**
```bash
python fem_crack_simulation_newton_raphson.py
```

**Purpose**: Simulates crack propagation driven by mechanical loading only
- **Input**: `robust_crack.npz` 
- **Workflow Documentation**: [`pure_mechanical.md`](pure_mechanical.md)
- **Algorithm**: Coupled mechanical-phase field system with penalty method
- **Loading**: Prescribed shear displacement at top boundary
- **Difference**: Direct method always converges; Newton-Raphson is exact nonlinear but may fail to converge
- **Output**: Time series visualization of displacement, phase field, and strain energy

### Step 2B: Electro-Mechanical Simulation  
```bash
python fem_electric_simulation.py
```

**Purpose**: Simulates crack propagation driven by combined electrical and mechanical forces
- **Input**: `robust_crack.msh`
- **Workflow Documentation**: [`electric-mechanical.md`](electric-mechanical.md) 
- **Algorithm**: 2-step solution (Electric → Mechanical+Phase field coupled)
- **Loading**: Current injection + mechanical displacement
- **Output**: Multi-field visualization (electric potential, displacement, energy density, phase field)

## 📁 Repository Structure

```
├── README.md                             # This file
├── robust_crack_mesh.py                 # Mesh generation system
├── fem_crack_simulation_updated.py      # Pure mechanical simulation (direct method)
├── fem_crack_simulation_newton_raphson.py # Pure mechanical simulation (Newton-Raphson)
├── fem_electric_simulation.py           # Electro-mechanical simulation  
├── pure_mechanical.md                   # Pure mechanical workflow documentation
├── electric-mechanical.md               # Electro-mechanical workflow documentation

## 🔧 Technical Specifications

### Material Properties
- **Steel-like material**: E = 210 GPa, ν = 0.3, ρ = 7380 kg/m³
- **Fracture energy**: Gc = 4000 N/m
- **Length scale**: l₀ = 12.5 μm
- **Body force**: Gravity (72.4×10³ N/m³ downward)

### Numerical Parameters
- **Element type**: 3-node linear triangles
- **Time step**: Δt = 1.0 s (quasi-static)
- **Integration**: 3-point Gauss quadrature
- **Solver**: Direct sparse (UMFPACK)
- **Solution Method**: Direct stiffness matrix (no iteration required)

### Boundary Conditions

#### Mechanical (Both Simulations)
- **Bottom edge**: Fixed displacement (u = 0)
- **Top edge**: Prescribed shear (uₓ = 1×10⁻⁵t mm, uᵧ = 0)
- **Left/Right edges**: Stress-free with minimal rigid body constraints

#### Electrical (Electro-Mechanical Only)
- **Left edge (above crack)**: Current injection (j₀ = 1×10⁸ A/m²)
- **Left edge (below crack)**: Current extraction (-j₀)
- **Right edge**: Ground (V = 0)
- **Top/Bottom edges**: Insulated (∂V/∂n = 0)

## 🚀 Quick Start

### Prerequisites
```bash
pip install numpy scipy matplotlib gmsh meshio
```

### Basic Usage

1. **Generate mesh**:
   ```bash
   python robust_crack_mesh.py
   ```
   This creates `robust_crack.msh` and `robust_crack.npz` files.

2. **Run mechanical simulation**:
   ```bash
   python fem_crack_simulation_updated.py
   ```
   Results saved in same directory with the excution code.

3. **Run electro-mechanical simulation**:
   ```bash
   python fem_electric_simulation.py
   # Select: 1 (test mode) or 2 (full simulation)
   ```
   Results saved as `electric_fem_step_*.png` files.

## 📊 Output and Visualization

### Pure Mechanical Simulation
- **Displacement magnitude** (μm)
- **Phase field** (crack damage indicator, φ ∈ [0,1])
- **Strain energy density** (J/m³, log scale)
- **Combined displacement + crack visualization**

### Electro-Mechanical Simulation  
- **Electric potential** (V)
- **Displacement magnitude** (μm)
- **Total energy density** (mechanical + electrical, J/m³)
- **Phase field evolution** (φ ∈ [0,1])

## 🧮 Mathematical Framework

### Phase Field Evolution
Both simulations solve the phase field equation with penalty-based irreversibility:

```
Gc l₀ ∇²φ - (Gc/l₀)φ + 2ψ(ε) + η/(m∆t)⟨φ̇⟩₋ᵐ = 2ψ(ε)
```

**Pure Mechanical**: ψ(ε) = elastic energy density  
**Electro-Mechanical**: ψ(ε) + W(E) = elastic + electric energy density

### Mechanical Equilibrium
Degraded elasticity with phase field coupling:
```
∇ · [g(φ) σ(ε)] + b = 0
g(φ) = (1-φ)² + k    (degradation function)
```

### Electric Potential (Electro-Mechanical)
Degraded conductivity equation:
```
∇ · [σ(φ) ∇V] = 0
σ(φ) = h(φ) σ₀    (degraded conductivity)
```

## ⚡ Performance Features

- **Optimized matrix assembly** with vectorized operations
- **Sparse matrix storage** (CSR format)
- **Efficient boundary condition enforcement**
- **Removed artificial stress concentrations** through proper constraint handling
- **Typical performance**: ~15 seconds per time step (5889 nodes, 17667 DOFs)

## 🔬 Key Features

### Mesh Quality
- **Adaptive refinement**: 5 μm elements at crack tip
- **Smooth transitions**: 10 μm → 25 μm element progression  
- **Crack void representation**: Direct geometric modeling (no artificial φ = 1 regions)
- **High aspect ratio handling**: Robust for thin crack geometry

### Physical Accuracy
- **Irreversibility constraint**: φ̇ ≥ 0 enforced via penalty method
- **Energy balance**: Conservative energy treatment
- **Stress-free crack surfaces**: Natural boundary conditions
- **Proper rigid body constraints**: Minimal artificial constraints

### Numerical Stability
- **Penalty method**: Robust irreversibility enforcement (η = 4×10⁶)
- **Residual stiffness**: Prevents singularity (k = 1×10⁻⁶)
- **Phase field bounds**: φ ∈ [0,1] enforced
- **Mesh-independent crack representation**

## 📚 Documentation

- **[`pure_mechanical.md`](pure_mechanical.md)**: Complete mathematical formulation for pure mechanical simulation
- **[`electric-mechanical.md`](electric-mechanical.md)**: Complete formulation for electro-mechanical coupling
- **Code comments**: Extensive inline documentation linking implementation to theory

## 🔧 Customization

### Simulation Parameters
Edit material properties, time steps, and simulation duration in the respective simulation files:
- `fem_crack_simulation_updated.py` (lines 36-56)
- `fem_electric_simulation.py` (lines 37-74)

### Mesh Parameters  
Modify crack geometry and mesh refinement in:
- `robust_crack_mesh.py` (lines 15-28)

### Loading Conditions
Adjust boundary conditions in the `apply_boundary_conditions()` methods of each simulation.

## 🐛 Troubleshooting

### Common Issues
1. **Missing mesh files**: Run `robust_crack_mesh.py` first
2. **Memory issues**: Reduce simulation time or mesh density
3. **Solution problems**: Check time step size and penalty parameters (no convergence issues with direct method)
4. **Performance**: Use test mode first (10 steps vs 2000+ steps)

### Performance Optimization
- **Matrix warnings**: Normal for LIL matrix operations during assembly
- **Step timing**: ~15s per step is expected for high-resolution mesh
- **Memory usage**: ~2-3 GB for typical simulation

## 📄 Citation

If you use this code in your research, please cite:

```bibtex
@software{phase_field_crack_simulation,
  title={crack-propagation},
  author={[philxie48]},
  year={2025},
  url={https://github.com/philxie48/crack-propagation}
}
```

## 📜 License

MIT license

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📞 Support

- **Issues**: Use GitHub Issues for bug reports and feature requests
- **Documentation**: Refer to the detailed workflow files (`pure_mechanical.md`,`electric-mechanical.md`)


---

**Keywords**: Phase field modeling, Crack propagation, Finite element method, Electro-mechanical coupling, Fracture mechanics, Computational materials science
