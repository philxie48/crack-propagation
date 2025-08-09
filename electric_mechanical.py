#!/usr/bin/env python3
"""
FEM Electro-Mechanical Crack Propagation Simulation
Based on fem-electric.md workflow

Solves coupled:
1. Electric potential (with degraded conductivity)
2. Mechanical displacement (with electrical energy effects) 
3. Phase field evolution (driven by mechanical + electrical energy)

Author: Generated from fem-electric.md workflow
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import time
import os
import meshio

class FEMElectricCrackSimulation:
    def __init__(self, mesh_file='robust_crack.msh', test_mode=False):
        """
        Initialize the electro-mechanical crack simulation
        
        Parameters:
        mesh_file: Path to mesh file (.msh format)
        test_mode: If True, run short test (10 steps), else full simulation (2000 steps)
        """
        self.mesh_file = mesh_file
        self.test_mode = test_mode
        
        # Load mesh
        self.load_mesh()
        
        # Material properties
        self.E = 210e9     # Pa (Young's modulus)
        self.nu = 0.3      # Poisson's ratio
        self.rho = 7.38e3  # kg/m³ (density)
        
        # Phase field properties
        self.l0 = 1.25e-5    # m (length scale) - standardized with mechanical simulation
        self.k = 1e-6      # Residual stiffness
        
        # Electrical properties
        self.sigma0 = 8.7e6  # S/m (electrical conductivity)
        self.eps0 = 8.854e-12  # F/m (vacuum permittivity)
        self.B = 1.24e7    # Constant B (dimensionless)
        self.Gcr_s = 4000  # N/m (critical energy release rate)
        self.h_bar = 1e-6  # Residual electrical conductivity
        
        # Penalty method parameters
        self.eta = 4e6     # Penalty parameter
        self.m = 2         # Penalty exponent
        
        # Current density (SI): 1e4 A/cm² = 1e8 A/m²
        self.j0 = 1e8      # A/m²
        
        # Body force (gravity)
        self.g = 9.81      # m/s²
        self.body_force = np.array([0, -72.4e3])  # N/m³ downward (from fem_em.md)
        
        # Time parameters
        if test_mode:
            self.dt = 1.0       # s (time step)
            self.t_final = 10.0 # s (final time - test)
            self.save_interval = 1  # Save every step
        else:
            self.dt = 1.0       # s (time step)
            self.t_final = 2000.0 # s (final time)
            self.save_interval = 20  # Save every 20 steps
            
        self.n_steps = int(self.t_final / self.dt)
        self.current_time = 0.0
        
        # Constitutive matrix (plane strain)
        self.D = self.compute_constitutive_matrix()
        
        # Initialize solution fields
        self.V = np.zeros(self.n_nodes)       # Electric potential
        self.u = np.zeros(2 * self.n_nodes)  # Displacement [ux1, uy1, ux2, uy2, ...]
        self.phi = np.zeros(self.n_nodes)    # Phase field
        self.phi_old = np.zeros(self.n_nodes)  # Previous time step phase field
        
        # Initialize boundary conditions
        self.find_boundary_nodes()
        
        print(f"Electro-mechanical simulation initialized:")
        print(f"  Mesh: {self.n_nodes} nodes, {self.n_elements} elements")
        print(f"  Mode: {'Test' if test_mode else 'Full'} simulation")
        print(f"  Time: {self.dt}s steps, {self.n_steps} total steps")
        print(f"  Current density: {self.j0:.1e} A/m²")
        
    def load_mesh(self):
        """Load mesh from .msh file"""
        mesh = meshio.read(self.mesh_file)
        
        # Extract coordinates and connectivity
        self.coordinates = mesh.points[:, :2]  # Only x, y coordinates
        self.n_nodes = len(self.coordinates)
        
        # Extract triangular elements (handle different meshio versions)
        triangles = None
        for cell_block in mesh.cells:
            if hasattr(cell_block, 'type') and hasattr(cell_block, 'data'):
                # New meshio version with CellBlock objects
                if cell_block.type == "triangle":
                    triangles = cell_block.data
                    break
            else:
                # Old meshio version with tuples
                cell_type, cell_data = cell_block
                if cell_type == "triangle":
                    triangles = cell_data
                    break
                
        if triangles is None:
            raise ValueError("No triangular elements found in mesh")
            
        self.elements = triangles
        self.n_elements = len(self.elements)
        
        print(f"Loaded mesh: {self.n_nodes} nodes, {self.n_elements} triangular elements")
        
    def find_boundary_nodes(self):
        """Find boundary nodes for different boundary conditions"""
        coords = self.coordinates
        tol = 1e-8
        
        # Find boundary nodes
        self.bottom_nodes = np.where(np.abs(coords[:, 1] - 0.0) < tol)[0]
        self.top_nodes = np.where(np.abs(coords[:, 1] - 1e-3) < tol)[0]  # 1mm
        self.left_nodes = np.where(np.abs(coords[:, 0] - 0.0) < tol)[0]
        self.right_nodes = np.where(np.abs(coords[:, 0] - 1e-3) < tol)[0]  # 1mm
        
        # Split left edge for electrical boundary conditions
        # Current out: y ∈ (0, 0.496mm), Current in: y ∈ (0.504mm, 1mm)
        left_coords = coords[self.left_nodes]
        self.current_out_nodes = self.left_nodes[
            (left_coords[:, 1] > 0.0) & (left_coords[:, 1] < 0.496e-3)
        ]
        self.current_in_nodes = self.left_nodes[
            (left_coords[:, 1] > 0.504e-3) & (left_coords[:, 1] < 1e-3)
        ]
        
        print(f"Boundary nodes found:")
        print(f"  Bottom: {len(self.bottom_nodes)}, Top: {len(self.top_nodes)}")
        print(f"  Left: {len(self.left_nodes)}, Right: {len(self.right_nodes)}")
        print(f"  Current out: {len(self.current_out_nodes)}, Current in: {len(self.current_in_nodes)}")
        
    def compute_constitutive_matrix(self):
        """Compute plane strain constitutive matrix (from fem_em.md)"""
        factor = self.E * (1 - self.nu) / ((1 + self.nu) * (1 - 2 * self.nu))
        D = np.array([
            [1.0, self.nu/(1-self.nu), 0],
            [self.nu/(1-self.nu), 1.0, 0],
            [0, 0, (1-2*self.nu)/(2*(1-self.nu))]
        ]) * factor
        return D
        
    def degradation_g(self, phi):
        """Mechanical degradation function g(φ) = (1-φ)² + k"""
        return (1 - phi)**2 + self.k
        
    def degradation_h(self, phi):
        """Electrical degradation function h(φ) = (1-φ)² + h̄"""
        return (1 - phi)**2 + self.h_bar
        
    def elastic_energy_density(self, strain):
        """Compute elastic energy density ψ(ε)"""
        if strain.ndim == 1:
            return 0.5 * strain @ self.D @ strain
        else:
            return 0.5 * np.sum((strain @ self.D) * strain, axis=1)
            
    def compute_B_matrices(self, element_nodes):
        """Compute B matrices for element (Jacobian-based method)"""
        # Element coordinates
        x1, y1 = element_nodes[0]
        x2, y2 = element_nodes[1]
        x3, y3 = element_nodes[2]
        
        # Jacobian matrix
        J = np.array([
            [x2 - x1, x3 - x1],
            [y2 - y1, y3 - y1]
        ])
        
        det_J = np.linalg.det(J)
        if abs(det_J) < 1e-12:
            print(f"Warning: Small Jacobian determinant: {det_J}")
            
        J_inv = np.linalg.inv(J)
        
        # Shape function derivatives in physical coordinates
        dN_dx = J_inv @ np.array([[-1, 1, 0], [-1, 0, 1]])
        
        # B matrix for displacement (strain-displacement)
        B_u = np.zeros((3, 6))  # 3 strains, 6 DOFs (2 per node)
        for i in range(3):
            B_u[0, 2*i] = dN_dx[0, i]      # ∂N_i/∂x for εxx
            B_u[1, 2*i+1] = dN_dx[1, i]    # ∂N_i/∂y for εyy
            B_u[2, 2*i] = dN_dx[1, i]      # ∂N_i/∂y for γxy
            B_u[2, 2*i+1] = dN_dx[0, i]    # ∂N_i/∂x for γxy
        
        # B matrix for gradients (electric field, phase field)
        B_grad = np.zeros((2, 3))  # 2 gradients, 3 nodes
        B_grad[0, :] = dN_dx[0, :]  # ∂N_i/∂x
        B_grad[1, :] = dN_dx[1, :]  # ∂N_i/∂y
        
        return B_u, B_grad, det_J
            
    def electric_energy_density(self, E_field):
        """Compute electric energy density W(E) = ½ε₀E²"""
        if E_field.ndim == 1:
            return 0.5 * self.eps0 * np.sum(E_field**2)
        else:
            return 0.5 * self.eps0 * np.sum(E_field**2, axis=1)
            
    def solve_electric_potential(self):
        """Step 1: Solve electric potential with degraded conductivity"""
        print("  Solving electric potential...")
        
        # Assemble conductivity matrix
        K_V = sp.lil_matrix((self.n_nodes, self.n_nodes))
        F_V = np.zeros(self.n_nodes)
        
        for elem_idx, element in enumerate(self.elements):
            # Element coordinates
            elem_coords = self.coordinates[element]
            
            # Use proper B-matrix computation
            B_u, B_grad, det_J = self.compute_B_matrices(elem_coords)
            area = 0.5 * abs(det_J)
            if area < 1e-16:
                continue
            
            # Element phase field values (for interpolation)
            phi_elem_values = self.phi[element]
            
            # Initialize element matrix
            K_elem = np.zeros((3, 3))
            
            # Gauss integration for conductivity matrix (using 3-point integration)
            gauss_points = [
                (1/6, 1/6, 1/3),      # Point 1
                (2/3, 1/6, 1/3),      # Point 2  
                (1/6, 2/3, 1/3)       # Point 3
            ]
            
            for xi, eta, weight in gauss_points:
                # Interpolate phase field at Gauss point
                N = np.array([1 - xi - eta, xi, eta])
                phi_gp = N.dot(phi_elem_values)
                
                # Degraded conductivity at Gauss point
                sigma_gp = self.sigma0 * self.degradation_h(phi_gp)
                
                # Add contribution to element matrix
                K_elem += area * weight * sigma_gp * (B_grad.T @ B_grad)
            
            # Assemble into global matrix
            for i in range(3):
                for j in range(3):
                    K_V[element[i], element[j]] += K_elem[i, j]
        
        # Convert to CSR for efficiency
        K_V = K_V.tocsr()
        
        # Apply electrical boundary conditions with length-weighted current distribution
        
        # Current injection at left edge (above crack) - Length-weighted distribution
        if len(self.current_in_nodes) > 0:
            in_coords = self.coordinates[self.current_in_nodes]
            if len(in_coords) > 1:
                # Sort nodes by y-coordinate
                sorted_idx = np.argsort(in_coords[:, 1])
                sorted_nodes = self.current_in_nodes[sorted_idx]
                sorted_coords = in_coords[sorted_idx]
                
                # Compute segment lengths between consecutive nodes
                segment_lengths = []
                for k in range(len(sorted_coords) - 1):
                    length = np.linalg.norm(sorted_coords[k+1] - sorted_coords[k])
                    segment_lengths.append(length)
                
                total_boundary_length = sum(segment_lengths)
                
                # Distribute current using length-weighted approach
                for i, node in enumerate(sorted_nodes):
                    if i == 0:
                        # First node: gets half of first segment
                        nodal_length = segment_lengths[0] / 2 if len(segment_lengths) > 0 else 0
                    elif i == len(sorted_nodes) - 1:
                        # Last node: gets half of last segment
                        nodal_length = segment_lengths[-1] / 2 if len(segment_lengths) > 0 else 0
                    else:
                        # Interior nodes: get half of adjacent segments
                        nodal_length = (segment_lengths[i-1] + segment_lengths[i]) / 2
                    
                    # Apply current density over the length represented by this node
                    current = self.j0 * nodal_length   # A/m² × m = A
                    F_V[node] += current
            else:
                # Single node: assume small area (1mm² = 1e-6 m²)
                area_m2 = 1e-6  # 1 mm² in m²
                current = self.j0 * area_m2  # A/m² × m² = A
                F_V[self.current_in_nodes[0]] += current
        
        # Current extraction at left edge (below crack) - Length-weighted distribution
        if len(self.current_out_nodes) > 0:
            out_coords = self.coordinates[self.current_out_nodes]
            if len(out_coords) > 1:
                # Sort nodes by y-coordinate
                sorted_idx = np.argsort(out_coords[:, 1])
                sorted_nodes = self.current_out_nodes[sorted_idx]
                sorted_coords = out_coords[sorted_idx]
                
                # Compute segment lengths between consecutive nodes
                segment_lengths = []
                for k in range(len(sorted_coords) - 1):
                    length = np.linalg.norm(sorted_coords[k+1] - sorted_coords[k])
                    segment_lengths.append(length)
                
                total_boundary_length = sum(segment_lengths)
                
                # Distribute current using length-weighted approach
                for i, node in enumerate(sorted_nodes):
                    if i == 0:
                        # First node: gets half of first segment
                        nodal_length = segment_lengths[0] / 2 if len(segment_lengths) > 0 else 0
                    elif i == len(sorted_nodes) - 1:
                        # Last node: gets half of last segment
                        nodal_length = segment_lengths[-1] / 2 if len(segment_lengths) > 0 else 0
                    else:
                        # Interior nodes: get half of adjacent segments
                        nodal_length = (segment_lengths[i-1] + segment_lengths[i]) / 2
                    
                    # Apply current density over the length represented by this node
                    current = self.j0 * nodal_length   # A/m² × m = A
                    F_V[node] -= current
            else:
                # Single node: assume small area (1mm² = 1e-6 m²)
                area_m2 = 1e-6  # 1 mm² in m²
                current = self.j0 * area_m2  # A/m² × m² = A
                F_V[self.current_out_nodes[0]] -= current
            
        # Ground at right edge: V = 0
        for node in self.right_nodes:
            # Set diagonal to 1, force to 0
            K_V[node, :] = 0
            K_V[node, node] = 1
            F_V[node] = 0
            
        # Solve for electric potential
        try:
            self.V = spsolve(K_V, F_V)
        except Exception as e:
            print(f"    Error solving electric potential: {e}")
            return False
            
        print(f"    Electric potential range: [{np.min(self.V):.6f}, {np.max(self.V):.6f}] V")
        return True
        

        

        
    def penalty_function(self, phi_dot):
        """Compute penalty function for irreversibility"""
        if phi_dot < 0:  # Healing attempt
            return self.eta / (self.m * self.dt) * (-phi_dot) ** self.m
        else:
            return 0
    
    def penalty_derivative(self, phi_dot):
        """Derivative of penalty function"""
        if phi_dot < 0:  # Healing attempt
            return self.eta / self.dt * (-phi_dot) ** (self.m - 1)
        else:
            return 0
            
    def shape_functions(self, xi, eta):
        """Linear triangle shape functions and derivatives"""
        N = np.array([1 - xi - eta, xi, eta])
        
        # Shape function derivatives in natural coordinates
        dN_dxi = np.array([-1, 1, 0])
        dN_deta = np.array([-1, 0, 1])
        
        return N, dN_dxi, dN_deta
        
    def solve_coupled_system(self):
        """Solve the three-field coupled system: Electric (sequential) + Mechanical-Phase (coupled)"""
        success = True
        # Step 1: Solve electric potential
        success &= self.solve_electric_potential()
        # Step 2: Solve mechanical displacement and phase field together (coupled)
        success &= self.solve_mechanical_phase_coupled()
        return success
        
    def solve_mechanical_phase_coupled(self):
        """Step 2: Solve mechanical displacement and phase field together (coupled with electric energy)"""
        print("  Solving coupled mechanical-phase field system...")
        
        # Electric energy will be computed directly at integration points - no pre-computation needed
        
        # Assemble coupled mechanical-phase field system
        n_dof_u = 2 * self.n_nodes  # Displacement DOFs
        n_dof_phi = self.n_nodes     # Phase field DOFs
        total_dof = n_dof_u + n_dof_phi
        
        # Initialize matrices
        K = sp.lil_matrix((total_dof, total_dof))
        F = np.zeros(total_dof)
        
        # Gauss points for triangle (3-point integration for better accuracy)
        gauss_points = [
            (1/6, 1/6, 1/3),      # Point 1
            (2/3, 1/6, 1/3),      # Point 2  
            (1/6, 2/3, 1/3)       # Point 3
        ]
        
        for elem_idx, element in enumerate(self.elements):
            elem_coords = self.coordinates[element]
            
            # Get B matrices
            B_u, B_grad, det_J = self.compute_B_matrices(elem_coords)
            area = 0.5 * abs(det_J)
            if area < 1e-16:
                continue
            
            # Element phase field values
            phi_elem = self.phi[element]
            phi_old_elem = self.phi_old[element]
            
            # Element displacement
            u_elem = np.zeros(6)
            for i in range(3):
                u_elem[2*i] = self.u[2*element[i]]      # ux
                u_elem[2*i+1] = self.u[2*element[i]+1]  # uy
            
            # EFFICIENCY FIX: Compute strain and energy density once per element
            # For linear triangles, strain is constant throughout the element
            strain = B_u @ u_elem
            psi_elem = self.elastic_energy_density(strain)  # Element-constant mechanical energy density
            
            # Gauss integration
            for xi, eta, weight in gauss_points:
                N, _, _ = self.shape_functions(xi, eta)
                
                # Interpolate phase field at Gauss point
                phi_gp = N.dot(phi_elem)  # φ at this integration point
                
                # Degradation function g(φ) = (1-φ)² + k (now spatially varying)
                g_phi = self.degradation_g(phi_gp)
                
                # Compute electric energy DIRECTLY at integration point
                V_elem = self.V[element]                                    # Nodal voltages
                E_field_gp = -B_grad @ V_elem                              # E-field at GP (element-constant)
                W_gp = 0.5 * self.eps0 * np.sum(E_field_gp**2)           # Energy at GP (direct)
                
                # Calculate voltage-dependent Gc at each integration point
                V_point = np.dot(N, self.V[element])  # Voltage at integration point
                Gc_effective = max(self.Gcr_s - self.B * V_point**2, 0.1 * self.Gcr_s)  # Pointwise Gc
                
                # K_uu: Mechanical stiffness matrix
                K_uu_elem = g_phi * B_u.T @ self.D @ B_u * area * weight
                
                # K_phi_phi: Phase field stiffness matrix (with electric energy coupling, now element-constant psi)
                # From fem-electric.md: Gc*l₀*∇φ·∇φ + (Gc/l₀ + 2ψ + 2W)*φ
                K_phi_phi_elem = (Gc_effective * self.l0 * B_grad.T @ B_grad + 
                                 (Gc_effective / self.l0 + 2 * psi_elem + 2 * W_gp) * np.outer(N, N)) * area * weight
                
                # K_u_phi: Coupling terms (mechanical effect on phase field)
                K_u_phi_elem = np.zeros((6, 3))
                dg_dphi = -2 * (1 - phi_gp)  # Derivative of degradation function (now pointwise)
                stress = self.D @ strain
                
                # K_phi_u: Coupling terms (phase field effect on mechanics)  
                K_phi_u_elem = np.zeros((3, 6))
                
                # Coupling terms (proper physics-based coupling, no arbitrary scaling, now element-constant psi)
                for i in range(3):
                    for j in range(3):
                        # Coupling between displacement and phase field nodes
                        K_u_phi_elem[2*i, j] += dg_dphi * stress[0] * N[j] * area * weight
                        K_u_phi_elem[2*i+1, j] += dg_dphi * stress[1] * N[j] * area * weight
                        
                        # Phase field effect on mechanics (include element-constant mechanical energy + pointwise electric energy)
                        K_phi_u_elem[i, 2*j] += 2 * (psi_elem + W_gp) * N[i] * area * weight
                        K_phi_u_elem[i, 2*j+1] += 2 * (psi_elem + W_gp) * N[i] * area * weight
                
                # Add penalty terms to K_phi_phi
                for i in range(3):
                    phi_dot = (phi_elem[i] - phi_old_elem[i]) / self.dt
                    penalty_deriv = self.penalty_derivative(phi_dot)
                    K_phi_phi_elem[i, i] += penalty_deriv * area * weight / self.dt
                
                # Assemble into global matrix with proper integer indexing
                # K_uu block
                for i in range(3):
                    for j in range(3):
                        for di in range(2):
                            for dj in range(2):
                                row = int(2 * element[i] + di)
                                col = int(2 * element[j] + dj)
                                K[row, col] += K_uu_elem[2*i + di, 2*j + dj]
                
                # K_phi_phi block
                for i in range(3):
                    for j in range(3):
                        row = int(n_dof_u + element[i])
                        col = int(n_dof_u + element[j])
                        K[row, col] += K_phi_phi_elem[i, j]
                
                # K_u_phi block
                for i in range(3):
                    for j in range(3):
                        for di in range(2):
                            row = int(2 * element[i] + di)
                            col = int(n_dof_u + element[j])
                            K[row, col] += K_u_phi_elem[2*i + di, j]
                
                # K_phi_u block
                for i in range(3):
                    for j in range(3):
                        for dj in range(2):
                            row = int(n_dof_u + element[i])
                            col = int(2 * element[j] + dj)
                            K[row, col] += K_phi_u_elem[i, 2*j + dj]
                
                # Force vector
                # Mechanical body force
                for i in range(3):
                    F[int(2*element[i])] += self.body_force[0] * N[i] * area * weight
                    F[int(2*element[i]+1)] += self.body_force[1] * N[i] * area * weight
                
                # Phase field driving force (with element-constant mechanical energy + pointwise electric energy and voltage-dependent Gc)
                for i in range(3):
                    # Force from energy: 2*(ψ + W) with element-constant psi + pointwise W
                    F[int(n_dof_u + element[i])] += 2 * (psi_elem + W_gp) * N[i] * area * weight
                    
                    # Penalty force (irreversibility constraint)
                    phi_dot = (phi_elem[i] - phi_old_elem[i]) / self.dt
                    if phi_dot < 0:  # Only apply penalty for healing attempts
                        penalty_force = self.eta / (self.m * self.dt) * (-phi_dot) ** self.m
                        F[int(n_dof_u + element[i])] -= penalty_force * N[i] * area * weight
        
        K = K.tocsr()
        
        # Apply boundary conditions
        K, F = self.apply_mechanical_phase_boundary_conditions(K, F)
        
        # Solve coupled system
        try:
            solution = spsolve(K, F)
        except Exception as e:
            print(f"    Error solving coupled system: {e}")
            return False
            
        # Extract solutions
        self.u = solution[:n_dof_u]
        phi_new = solution[n_dof_u:]
        
        # Update phase field with bounds and irreversibility
        self.phi_old = self.phi.copy()
        phi_new = np.clip(phi_new, 0, 1)  # Phase field bounds
        
        # Enforce irreversibility: φ can only increase
        self.phi = np.maximum(self.phi, phi_new)
        
        # Statistics
        max_u = np.max(np.abs(self.u))
        max_phi = np.max(self.phi)
        min_phi = np.min(self.phi)
        crack_area = np.sum(self.phi > 0.1) / self.n_nodes
        max_V = np.max(np.abs(self.V))
        
        # Compute effective Gc statistics
        V_mean = np.mean(np.abs(self.V))
        Gc_min = self.Gcr_s - self.B * np.max(self.V)**2
        Gc_reduction = (self.Gcr_s - Gc_min) / self.Gcr_s * 100
        
        # Energy statistics - compute max electric energy over elements
        max_electric_energy = 0.0
        for elem_idx, element in enumerate(self.elements):
            elem_coords = self.coordinates[element]
            B_u, B_grad, det_J = self.compute_B_matrices(elem_coords)
            if abs(det_J) < 1e-16:
                continue
            V_elem = self.V[element]
            E_field = -B_grad @ V_elem
            W_elem = 0.5 * self.eps0 * np.sum(E_field**2)
            max_electric_energy = max(max_electric_energy, W_elem)
        
        print(f"    Max displacement: {max_u*1e6:.3f} μm")
        print(f"    Phase field range: [{min_phi:.6f}, {max_phi:.6f}]")
        print(f"    Crack area fraction: {crack_area:.4f}")
        print(f"    Max voltage: {max_V:.6f} V, Mean |V|: {V_mean:.6f} V")
        print(f"    Gc reduction: {Gc_reduction:.2f}% (min Gc = {max(Gc_min, 0.1*self.Gcr_s):.0f} N/m)")
        print(f"    Electric energy contribution: {max_electric_energy:.2e} J/m³")
        
        return True
        
    def apply_mechanical_phase_boundary_conditions(self, K, F):
        """Apply boundary conditions for coupled mechanical-phase field system"""
        n_dof_u = 2 * self.n_nodes
        
        # Bottom boundary: u = 0 (fixed)
        for node in self.bottom_nodes:
            # ux = 0
            row = int(2 * node)
            K[row, :] = 0
            K[row, row] = 1
            F[row] = 0
            
            # uy = 0
            row = int(2 * node + 1)
            K[row, :] = 0
            K[row, row] = 1
            F[row] = 0
        
        # Top boundary: prescribed shear ux = 1e-5 * t mm, uy = 0
        prescribed_ux = 1e-5 * self.current_time * 1e-3  # Convert mm to m
        
        for node in self.top_nodes:
            # ux = prescribed displacement
            row = int(2 * node)
            K[row, :] = 0
            K[row, row] = 1
            F[row] = prescribed_ux
            
            # uy = 0
            row = int(2 * node + 1)
            K[row, :] = 0
            K[row, row] = 1
            F[row] = 0
        
        # Add rigid body constraints to prevent spurious motion
        # Constraint 1: Fix one point completely (already done with bottom boundary)
        
        # Constraint 2: Prevent rotation - fix x-displacement of one left boundary node
        if len(self.left_nodes) > 0:
            constraint_node = self.left_nodes[len(self.left_nodes)//2]  # Middle left node
            if constraint_node not in self.bottom_nodes:  # Avoid double constraint
                row = int(2 * constraint_node)  # ux = 0
                K[row, :] = 0
                K[row, row] = 1
                F[row] = 0
        
        # Constraint 3: Prevent rigid body translation in y - constrain average y-displacement
        # of right boundary nodes (except those already constrained)
        free_right_nodes = [n for n in self.right_nodes if n not in self.bottom_nodes and n not in self.top_nodes]
        if len(free_right_nodes) > 1:
            # Add constraint: average uy of free right nodes = 0
            # For simplicity, just fix the middle right node's y-displacement
            constraint_node = free_right_nodes[len(free_right_nodes)//2]
            row = int(2 * constraint_node + 1)  # uy = 0
            K[row, :] = 0
            K[row, row] = 1
            F[row] = 0
        
        return K, F
        
    def visualize_results(self, step):
        """Visualize all four fields: V, |u|, ψ+W, φ"""
        if step % self.save_interval != 0 and step != self.n_steps:
            return
            
        print("  Creating visualization...")
        
        # Compute displacement magnitude
        u_mag = np.sqrt(self.u[::2]**2 + self.u[1::2]**2)
        
        # Compute strain energy density
        strain_energy = np.zeros(self.n_nodes)
        electric_energy = np.zeros(self.n_nodes)
        
        for elem_idx, element in enumerate(self.elements):
            elem_coords = self.coordinates[element]
            
            # Use proper B-matrix computation
            B_u, B_grad, det_J = self.compute_B_matrices(elem_coords)
            area = 0.5 * abs(det_J)
            if area < 1e-16:
                continue
                
            # Strain energy
            u_elem = np.zeros(6)
            for i in range(3):
                u_elem[2*i] = self.u[2*element[i]]
                u_elem[2*i+1] = self.u[2*element[i]+1]
            strain = B_u @ u_elem
            psi_elem = self.elastic_energy_density(strain)
            
            # Electric energy
            V_elem = self.V[element]
            E_field = -B_grad @ V_elem
            W_elem = self.electric_energy_density(E_field)
            
            # Distribute to nodes
            for i in range(3):
                strain_energy[element[i]] += psi_elem / 3
                electric_energy[element[i]] += W_elem / 3
        
        total_energy = strain_energy + electric_energy
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        coords_mm = self.coordinates * 1000  # Convert to mm
        
        # Plot 1: Electric potential
        im1 = axes[0,0].tripcolor(coords_mm[:, 0], coords_mm[:, 1], self.elements, 
                                  self.V, shading='flat', cmap='RdBu_r')
        axes[0,0].set_title(f'Electric Potential (V) - t={self.current_time:.0f}s')
        axes[0,0].set_xlabel('x (mm)')
        axes[0,0].set_ylabel('y (mm)')
        axes[0,0].set_aspect('equal')
        plt.colorbar(im1, ax=axes[0,0], label='Potential (V)')
        
        # Plot 2: Displacement magnitude
        im2 = axes[0,1].tripcolor(coords_mm[:, 0], coords_mm[:, 1], self.elements, 
                                  u_mag*1e6, shading='flat', cmap='viridis')
        axes[0,1].set_title(f'Displacement Magnitude - t={self.current_time:.0f}s')
        axes[0,1].set_xlabel('x (mm)')
        axes[0,1].set_ylabel('y (mm)')
        axes[0,1].set_aspect('equal')
        plt.colorbar(im2, ax=axes[0,1], label='|u| (μm)')
        
        # Plot 3: Total energy density (mechanical + electrical)
        im3 = axes[1,0].tripcolor(coords_mm[:, 0], coords_mm[:, 1], self.elements, 
                                  total_energy, shading='flat', cmap='plasma')
        axes[1,0].set_title(f'Total Energy Density (ψ+W) - t={self.current_time:.0f}s')
        axes[1,0].set_xlabel('x (mm)')
        axes[1,0].set_ylabel('y (mm)')
        axes[1,0].set_aspect('equal')
        plt.colorbar(im3, ax=axes[1,0], label='Energy (J/m³)')
        
        # Plot 4: Phase field
        im4 = axes[1,1].tripcolor(coords_mm[:, 0], coords_mm[:, 1], self.elements, 
                                  self.phi, shading='flat', cmap='Reds', vmin=0, vmax=1)
        axes[1,1].set_title(f'Phase Field (φ) - t={self.current_time:.0f}s')
        axes[1,1].set_xlabel('x (mm)')
        axes[1,1].set_ylabel('y (mm)')
        axes[1,1].set_aspect('equal')
        plt.colorbar(im4, ax=axes[1,1], label='φ')
        
        plt.tight_layout()
        filename = f'electric_fem_step_{step:04d}_t{self.current_time:04.0f}s.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        
        # Print statistics
        print(f"    Visualization saved: {filename}")
        print(f"    Max potential: {np.max(self.V):.3f} V")
        print(f"    Max displacement: {np.max(u_mag)*1e6:.3f} μm")
        print(f"    Max energy density: {np.max(total_energy):.2e} J/m³")
        print(f"    Max phase field: {np.max(self.phi):.6f}")
        
    def run_simulation(self):
        """Run the complete electro-mechanical simulation"""
        print("="*80)
        print("STARTING ELECTRO-MECHANICAL CRACK PROPAGATION SIMULATION")
        print("="*80)
        print(f"Mode: {'Test' if self.test_mode else 'Full'}")
        print(f"Time: {self.dt}s steps × {self.n_steps} = {self.t_final}s total")
        print(f"Current density: {self.j0:.1e} A/m²")
        print(f"Saving every {self.save_interval} steps")
        print("")
        
        start_time = time.time()
        
        for step in range(self.n_steps + 1):
            self.current_time = step * self.dt
            
            # Progress reporting
            if step % max(1, self.n_steps//20) == 0 or step % self.save_interval == 0:
                progress = step / self.n_steps * 100
                print(f"Step {step}/{self.n_steps} (t = {self.current_time:.0f}s) - {progress:.1f}% complete")
            
            # Solve coupled system
            step_start = time.time()
            success = self.solve_coupled_system()
            step_time = time.time() - step_start
            
            if not success:
                print(f"  Simulation failed at step {step}")
                break
                
            if step % self.save_interval == 0 or step == self.n_steps:
                print(f"  Step completed in {step_time:.3f}s")
            
            # Visualize results
            self.visualize_results(step)
            
            if step % max(1, self.n_steps//20) == 0:
                print("")
        
        total_time = time.time() - start_time
        print("="*80)
        print("ELECTRO-MECHANICAL SIMULATION COMPLETED")
        print("="*80)
        print(f"Total simulation time: {total_time:.2f}s")
        print(f"Average time per step: {total_time/self.n_steps:.3f}s")

def main():
    """Main function"""
    
    # Check if mesh file exists
    mesh_file = 'robust_crack.msh'
    if not os.path.exists(mesh_file):
        print(f"Error: Mesh file '{mesh_file}' not found!")
        print("Please ensure the robust_crack.msh file is in the current directory.")
        return
    
    # Ask user for simulation mode
    print("Electro-Mechanical FEM Crack Simulation")
    print("1. Test mode: 10 steps, save every step")
    print("2. Full mode: 2000 steps, save every 20 steps")
    
    mode = input("Select mode (1 or 2): ").strip()
    test_mode = (mode == "1")
    
    if test_mode:
        print("Running test simulation...")
    else:
        print("Running full simulation...")
    
    # Create and run simulation
    sim = FEMElectricCrackSimulation(mesh_file, test_mode)
    sim.run_simulation()

if __name__ == "__main__":
    main()