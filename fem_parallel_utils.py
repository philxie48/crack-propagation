#!/usr/bin/env python3
"""
Parallel Computing Utilities for FEM Crack Simulation
- Provides multi-core support for element assembly
- Automatic CPU detection and load balancing
- Memory-efficient parallel processing
"""

import numpy as np
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
import os
from functools import partial

class ParallelConfig:
    """Configuration for parallel computing"""
    
    def __init__(self, n_cores=None, verbose=True):
        """
        Initialize parallel configuration
        
        Parameters:
        -----------
        n_cores : int, optional
            Number of cores to use. If None, auto-detect optimal number
        verbose : bool
            Print configuration information
        """
        self.total_cores = cpu_count()
        
        if n_cores is None:
            # Auto-detect optimal number of cores
            self.n_cores = self._get_optimal_cores()
        else:
            self.n_cores = min(n_cores, self.total_cores)
        
        if verbose:
            print("="*70)
            print("PARALLEL COMPUTING CONFIGURATION")
            print("="*70)
            print(f"Total CPU cores available: {self.total_cores}")
            print(f"Cores to be used: {self.n_cores}")
            print(f"Parallel efficiency: {self.n_cores/self.total_cores*100:.1f}%")
            print(f"Recommended for: {'Small problems' if self.n_cores <= 4 else 'Medium problems' if self.n_cores <= 10 else 'Large problems'}")
            print("="*70)
    
    def _get_optimal_cores(self):
        """
        Determine optimal number of cores based on system
        
        Strategy:
        - Leave 1-2 cores free for system operations
        - Use 75-90% of available cores for computation
        """
        if self.total_cores <= 4:
            return max(2, self.total_cores - 1)
        elif self.total_cores <= 8:
            return self.total_cores - 2
        elif self.total_cores <= 16:
            return int(self.total_cores * 0.85)
        else:
            return int(self.total_cores * 0.90)
    
    def get_chunk_size(self, n_elements):
        """
        Calculate optimal chunk size for element processing
        
        Parameters:
        -----------
        n_elements : int
            Total number of elements
            
        Returns:
        --------
        chunk_size : int
            Number of elements per chunk
        """
        # Aim for at least 2-3 chunks per core for load balancing
        min_chunks = self.n_cores * 3
        chunk_size = max(1, n_elements // min_chunks)
        return chunk_size


def compute_element_residuals_batch(element_indices, nodes, elements, u, phi, phi_old, 
                                   D, Gc, l0, k, eta, m, dt, body_force):
    """
    Compute residuals for a batch of elements (worker function)
    
    This function is designed to be called in parallel by multiple processes.
    Each process handles a subset of elements.
    
    Parameters:
    -----------
    element_indices : array
        Indices of elements to process in this batch
    nodes : array
        Node coordinates
    elements : array
        Element connectivity
    u : array
        Displacement field
    phi : array
        Phase field
    phi_old : array
        Previous phase field
    D : array
        Constitutive matrix
    Gc, l0, k, eta, m, dt : float
        Material and numerical parameters
    body_force : array
        Body force vector
        
    Returns:
    --------
    R_u_batch : dict
        Dictionary mapping DOF indices to residual contributions
    R_phi_batch : dict
        Dictionary mapping node indices to residual contributions
    """
    # Fix for Windows multiprocessing + matplotlib issue
    import os
    os.environ['MPLBACKEND'] = 'Agg'  # Use non-interactive backend
    
    from scipy.linalg import inv, det
    import numpy as np
    
    # Initialize batch residuals as dictionaries for sparse assembly
    R_u_batch = {}
    R_phi_batch = {}
    
    # Gauss points for triangle
    gauss_points = [
        (1/6, 1/6, 1/3),
        (2/3, 1/6, 1/3),
        (1/6, 2/3, 1/3)
    ]
    
    for elem_id in element_indices:
        element = elements[elem_id]
        element_nodes = nodes[element]
        
        # Compute B matrices
        x1, y1 = element_nodes[0]
        x2, y2 = element_nodes[1]
        x3, y3 = element_nodes[2]
        
        J = np.array([[x2 - x1, x3 - x1],
                      [y2 - y1, y3 - y1]])
        det_J = det(J)
        area = 0.5 * abs(det_J)
        
        if abs(det_J) < 1e-12:
            continue
        
        J_inv = inv(J)
        dN_dx = J_inv @ np.array([[-1, 1, 0], [-1, 0, 1]])
        
        # B matrix for displacement
        B_u = np.zeros((3, 6))
        for i in range(3):
            B_u[0, 2*i] = dN_dx[0, i]
            B_u[1, 2*i+1] = dN_dx[1, i]
            B_u[2, 2*i] = dN_dx[1, i]
            B_u[2, 2*i+1] = dN_dx[0, i]
        
        # B matrix for phase field
        B_phi = np.zeros((2, 3))
        B_phi[0, :] = dN_dx[0, :]
        B_phi[1, :] = dN_dx[1, :]
        
        # Element phase field values
        phi_elem = phi[element]
        phi_old_elem = phi_old[element]
        
        # Element displacement
        u_elem = np.zeros(6)
        for i in range(3):
            u_elem[2*i] = u[int(2*element[i])]
            u_elem[2*i+1] = u[int(2*element[i]+1)]
        
        # Compute strain and stress
        strain = B_u @ u_elem
        stress = D @ strain
        psi_elem = max(0, 0.5 * strain @ stress)
        
        # Gauss integration
        for xi, eta, weight in gauss_points:
            N = np.array([1 - xi - eta, xi, eta])
            
            phi_gp = N.dot(phi_elem)
            g_phi = (1 - phi_gp)**2 + k
            
            # Mechanical residual
            R_u_internal = g_phi * B_u.T @ stress * area * weight
            R_u_body = np.zeros(6)
            for i in range(3):
                R_u_body[2*i] = body_force[0] * N[i] * area * weight
                R_u_body[2*i+1] = body_force[1] * N[i] * area * weight
            
            R_u_elem = R_u_internal - R_u_body
            
            # Phase field residual
            phi_grad = B_phi @ phi_elem
            phi_dot_elem = (phi_elem - phi_old_elem) / dt
            
            R_phi_diffusion = Gc * l0 * B_phi.T @ phi_grad * area * weight
            phi_at_gauss = N @ phi_elem
            R_phi_resistance = (Gc / l0 + 2 * psi_elem) * N * phi_at_gauss * area * weight
            R_phi_driving = -2 * psi_elem * N * area * weight
            
            R_phi_penalty = np.zeros(3)
            for i in range(3):
                phi_dot_i = phi_dot_elem[i]
                if phi_dot_i < 0:
                    penalty_value = eta / (m * dt) * (-phi_dot_i) ** m
                    R_phi_penalty[i] = -penalty_value * N[i] * area * weight
            
            R_phi_elem = R_phi_diffusion + R_phi_resistance + R_phi_driving + R_phi_penalty
            
            # Accumulate into batch dictionaries
            u_dofs = [2*element[0], 2*element[0]+1, 
                     2*element[1], 2*element[1]+1, 
                     2*element[2], 2*element[2]+1]
            
            for i, dof in enumerate(u_dofs):
                if dof not in R_u_batch:
                    R_u_batch[dof] = 0.0
                R_u_batch[dof] += R_u_elem[i]
            
            for i, node in enumerate(element):
                if node not in R_phi_batch:
                    R_phi_batch[node] = 0.0
                R_phi_batch[node] += R_phi_elem[i]
    
    return R_u_batch, R_phi_batch


def compute_residuals_parallel(nodes, elements, u, phi, phi_old, D, Gc, l0, k, 
                               eta, m, dt, body_force, n_cores=None):
    """
    Compute residuals using parallel processing
    
    Parameters:
    -----------
    nodes : array
        Node coordinates
    elements : array
        Element connectivity
    u : array
        Displacement field
    phi : array
        Phase field
    phi_old : array
        Previous phase field
    D : array
        Constitutive matrix
    Gc, l0, k, eta, m, dt : float
        Material and numerical parameters
    body_force : array
        Body force vector
    n_cores : int, optional
        Number of cores to use
        
    Returns:
    --------
    R_u : array
        Mechanical residual
    R_phi : array
        Phase field residual
    """
    config = ParallelConfig(n_cores=n_cores, verbose=False)
    n_elements = len(elements)
    n_nodes = len(nodes)
    
    # Split elements into chunks for parallel processing
    chunk_size = config.get_chunk_size(n_elements)
    element_chunks = [np.arange(i, min(i + chunk_size, n_elements)) 
                     for i in range(0, n_elements, chunk_size)]
    
    # Create worker function with fixed parameters
    worker = partial(compute_element_residuals_batch,
                    nodes=nodes, elements=elements, u=u, phi=phi, phi_old=phi_old,
                    D=D, Gc=Gc, l0=l0, k=k, eta=eta, m=m, dt=dt, body_force=body_force)
    
    # Process chunks in parallel
    with Pool(processes=config.n_cores) as pool:
        results = pool.map(worker, element_chunks)
    
    # Assemble results from all workers
    R_u = np.zeros(2 * n_nodes)
    R_phi = np.zeros(n_nodes)
    
    for R_u_batch, R_phi_batch in results:
        for dof, value in R_u_batch.items():
            R_u[int(dof)] += value
        for node, value in R_phi_batch.items():
            R_phi[int(node)] += value
    
    return R_u, R_phi


def benchmark_parallel_performance(n_elements_list=[100, 500, 1000, 5000], 
                                  n_cores_list=None):
    """
    Benchmark parallel performance for different problem sizes and core counts
    
    Parameters:
    -----------
    n_elements_list : list
        List of element counts to test
    n_cores_list : list, optional
        List of core counts to test. If None, test [1, 2, 4, 8, ...]
        
    Returns:
    --------
    results : dict
        Benchmark results
    """
    import time
    
    if n_cores_list is None:
        max_cores = cpu_count()
        n_cores_list = [1, 2, 4, 8, 16, max_cores]
        n_cores_list = [n for n in n_cores_list if n <= max_cores]
    
    print("="*70)
    print("PARALLEL PERFORMANCE BENCHMARK")
    print("="*70)
    print(f"Testing with {cpu_count()} available cores")
    print(f"Element counts: {n_elements_list}")
    print(f"Core counts: {n_cores_list}")
    print("="*70)
    
    results = {}
    
    for n_elements in n_elements_list:
        print(f"\nTesting with {n_elements} elements:")
        results[n_elements] = {}
        
        # Create dummy data
        nodes = np.random.rand(n_elements, 2)
        elements = np.random.randint(0, n_elements, (n_elements, 3))
        u = np.random.rand(2 * n_elements)
        phi = np.random.rand(n_elements)
        phi_old = np.random.rand(n_elements)
        D = np.eye(3)
        
        baseline_time = None
        
        for n_cores in n_cores_list:
            start = time.time()
            try:
                R_u, R_phi = compute_residuals_parallel(
                    nodes, elements, u, phi, phi_old, D,
                    Gc=4000, l0=1.25e-5, k=1e-6, eta=4e6, m=2, dt=1.0,
                    body_force=np.array([0, -72.4e3]),
                    n_cores=n_cores
                )
                elapsed = time.time() - start
                
                if baseline_time is None:
                    baseline_time = elapsed
                
                speedup = baseline_time / elapsed
                efficiency = speedup / n_cores * 100
                
                results[n_elements][n_cores] = {
                    'time': elapsed,
                    'speedup': speedup,
                    'efficiency': efficiency
                }
                
                print(f"  {n_cores:2d} cores: {elapsed:.4f}s | "
                      f"Speedup: {speedup:.2f}x | "
                      f"Efficiency: {efficiency:.1f}%")
            except Exception as e:
                print(f"  {n_cores:2d} cores: Failed - {e}")
    
    print("\n" + "="*70)
    print("RECOMMENDATIONS:")
    print("="*70)
    
    # Find optimal core count for largest problem
    if n_elements_list:
        largest = max(n_elements_list)
        if largest in results:
            best_cores = max(results[largest].keys(), 
                           key=lambda k: results[largest][k]['speedup'])
            best_speedup = results[largest][best_cores]['speedup']
            print(f"For {largest} elements:")
            print(f"  Optimal cores: {best_cores}")
            print(f"  Best speedup: {best_speedup:.2f}x")
            print(f"  Recommended for production: {best_cores} cores")
    
    return results


if __name__ == "__main__":
    # Run configuration check
    config = ParallelConfig()
    
    print("\nRunning quick benchmark...")
    benchmark_parallel_performance(n_elements_list=[500, 2000], 
                                  n_cores_list=[1, 2, 4, 8])

