#!/usr/bin/env python3
"""
Robust crack mesh generator with proper scaling
Uses appropriate mesh sizes and robust GMSH settings
"""

import gmsh
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class RobustCrackMesh:
    def __init__(self):
        # Domain parameters (in meters)
        self.Lx = 1.0e-3  # 1.0 mm
        self.Ly = 1.0e-3  # 1.0 mm
        
        # Crack parameters (in meters) - use reasonable ratios
        self.crack_length = 0.5e-3  # 0.5 mm
        self.crack_width = 8.0e-6  # 20.0 μm (larger for stable meshing)
        self.crack_center_y = 0.5e-3  # Center of domain
        
        # Mesh parameters (in meters) - properly scaled
        self.crack_tip_size = 5.0e-6   # Fine mesh at crack tip: 5 μm
        self.transition_size = 10.0e-6  # Transition size: 10 μm
        self.bulk_size = 25.0e-6       # Bulk mesh size: 25 μm
        self.tip_refinement_radius = 0.1e-3  # Refinement region radius: 0.1 mm
        
        # Crack geometry
        self.crack_tip_x = self.crack_length
        self.crack_tip_y = self.crack_center_y
        self.crack_y_bottom = self.crack_center_y - self.crack_width/2
        self.crack_y_top = self.crack_center_y + self.crack_width/2

    def create_robust_mesh(self, filename="robust_crack.msh"):
        """Create robust mesh using GMSH with proper geometry"""
        
        # Initialize GMSH
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("RobustCrackModel")
        
        print(f"Creating robust crack mesh:")
        print(f"  Domain: {self.Lx*1e3:.1f} × {self.Ly*1e3:.1f} mm")
        print(f"  Crack: from (0, {self.crack_center_y*1e3:.3f}) to ({self.crack_tip_x*1e3:.3f}, {self.crack_center_y*1e3:.3f}) mm")
        print(f"  Crack width: {self.crack_width*1e6:.1f} μm")
        print(f"  Mesh sizes: tip={self.crack_tip_size*1e6:.1f}μm, transition={self.transition_size*1e6:.1f}μm, bulk={self.bulk_size*1e6:.1f}μm")
        
        # Create outer domain rectangle
        rect = gmsh.model.occ.addRectangle(0, 0, 0, self.Lx, self.Ly)
        
        # Create crack void rectangle
        crack_rect = gmsh.model.occ.addRectangle(0, self.crack_y_bottom, 0, 
                                               self.crack_tip_x, self.crack_width)
        
        # Create domain with crack void (boolean difference)
        domain = gmsh.model.occ.cut([(2, rect)], [(2, crack_rect)], removeObject=True, removeTool=True)
        
        # Synchronize the CAD representation
        gmsh.model.occ.synchronize()
        
        # Get the entities after boolean operation
        domain_surfaces = gmsh.model.getEntities(2)
        domain_curves = gmsh.model.getEntities(1)
        domain_points = gmsh.model.getEntities(0)
        
        print(f"Created geometry with {len(domain_surfaces)} surfaces, {len(domain_curves)} curves, {len(domain_points)} points")
        
        # Set mesh sizes at key points
        # Find crack tip points
        all_points = gmsh.model.getEntities(0)
        for dim, tag in all_points:
            x, y, z = gmsh.model.getValue(0, tag, [])
            
            # Crack tip points (right edge of crack)
            if abs(x - self.crack_tip_x) < 1e-8 and abs(y - self.crack_center_y) < self.crack_width:
                gmsh.model.mesh.setSize([(0, tag)], self.crack_tip_size)
            # Crack left edge points
            elif abs(x) < 1e-8 and abs(y - self.crack_center_y) < self.crack_width:
                gmsh.model.mesh.setSize([(0, tag)], self.transition_size)
            # Domain corner points
            elif (abs(x) < 1e-8 or abs(x - self.Lx) < 1e-8) and (abs(y) < 1e-8 or abs(y - self.Ly) < 1e-8):
                gmsh.model.mesh.setSize([(0, tag)], self.bulk_size)
        
        # Create distance field from crack tip for refinement
        # Find crack tip curves
        crack_tip_curves = []
        for dim, tag in domain_curves:
            # Get curve bounds
            curve_points = gmsh.model.getBoundary([(1, tag)], combined=False, oriented=False)
            
            # Check if this curve is near the crack tip
            is_crack_tip_curve = False
            for pdim, ptag in curve_points:
                x, y, z = gmsh.model.getValue(0, ptag, [])
                if abs(x - self.crack_tip_x) < 1e-8 and abs(y - self.crack_center_y) < self.crack_width:
                    is_crack_tip_curve = True
                    break
            
            if is_crack_tip_curve:
                crack_tip_curves.append(tag)
        
        if crack_tip_curves:
            # Distance field from crack tip
            gmsh.model.mesh.field.add("Distance", 1)
            gmsh.model.mesh.field.setNumbers(1, "CurvesList", crack_tip_curves)
            
            # Threshold field for refinement
            gmsh.model.mesh.field.add("Threshold", 2)
            gmsh.model.mesh.field.setNumber(2, "InField", 1)
            gmsh.model.mesh.field.setNumber(2, "SizeMin", self.crack_tip_size)
            gmsh.model.mesh.field.setNumber(2, "SizeMax", self.bulk_size)
            gmsh.model.mesh.field.setNumber(2, "DistMin", 0)
            gmsh.model.mesh.field.setNumber(2, "DistMax", self.tip_refinement_radius)
            
            # Set as background mesh
            gmsh.model.mesh.field.setAsBackgroundMesh(2)
        
        # Set robust mesh options
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
        gmsh.option.setNumber("Mesh.RecombineAll", 0)  # Keep triangular elements
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.option.setNumber("Mesh.Smoothing", 3)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)  # Use point sizes
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.0)
        
        # Generate mesh
        print("\nGenerating mesh...")
        try:
            gmsh.model.mesh.generate(2)
            print("✓ 2D mesh generated successfully")
        except Exception as e:
            print(f"Error during mesh generation: {e}")
            gmsh.finalize()
            raise
        
        # Add physical groups
        gmsh.model.addPhysicalGroup(2, [s[1] for s in domain_surfaces], 1)
        gmsh.model.setPhysicalName(2, 1, "Material")
        
        # Get boundary curves for physical groups
        boundary_curves = []
        crack_curves = []
        
        for dim, tag in domain_curves:
            # Get curve bounds to classify
            curve_points = gmsh.model.getBoundary([(1, tag)], combined=False, oriented=False)
            
            # Check if this is a crack boundary or domain boundary
            is_crack_boundary = False
            for pdim, ptag in curve_points:
                x, y, z = gmsh.model.getValue(0, ptag, [])
                # If point is on crack boundary
                if (abs(y - self.crack_y_bottom) < 1e-8 or abs(y - self.crack_y_top) < 1e-8) and x <= self.crack_tip_x:
                    is_crack_boundary = True
                    break
                elif abs(x - self.crack_tip_x) < 1e-8 and self.crack_y_bottom <= y <= self.crack_y_top:
                    is_crack_boundary = True
                    break
            
            if is_crack_boundary:
                crack_curves.append(tag)
            else:
                boundary_curves.append(tag)
        
        if boundary_curves:
            gmsh.model.addPhysicalGroup(1, boundary_curves, 10)
            gmsh.model.setPhysicalName(1, 10, "DomainBoundary")
        
        if crack_curves:
            gmsh.model.addPhysicalGroup(1, crack_curves, 20)
            gmsh.model.setPhysicalName(1, 20, "CrackBoundary")
        
        # Save mesh
        gmsh.write(filename)
        print(f"✓ Mesh saved as {filename}")
        
        # Get mesh data
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements()
        
        # Extract triangle elements
        triangles = None
        for i, elem_type in enumerate(element_types):
            if elem_type == 2:  # Triangle
                triangles = element_node_tags[i].reshape(-1, 3) - 1
                break
        
        if triangles is None:
            triangles = np.array([])
        
        nodes = node_coords.reshape(-1, 3)[:, :2]
        
        print(f"Final mesh: {len(nodes)} nodes, {len(triangles)} triangular elements")
        
        gmsh.finalize()
        
        return nodes, triangles
    
    def plot_mesh(self, nodes, triangles, filename="robust_crack_mesh.png"):
        """Plot the robust mesh"""
        if len(triangles) == 0:
            print("No triangles to plot!")
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Plot 1: Full domain
        for triangle in triangles:
            triangle_nodes = nodes[triangle]
            triangle_closed = np.vstack([triangle_nodes, triangle_nodes[0]])
            ax1.plot(triangle_closed[:, 0]*1e3, triangle_closed[:, 1]*1e3, 'b-', linewidth=0.3, alpha=0.7)
        
        # Draw crack void
        crack_rect = patches.Rectangle((0, self.crack_y_bottom*1e3), 
                                     self.crack_tip_x*1e3, 
                                     self.crack_width*1e3,
                                     linewidth=2, edgecolor='red', 
                                     facecolor='white', alpha=0.9,
                                     label='Crack Void')
        ax1.add_patch(crack_rect)
        
        # Crack tip
        ax1.plot(self.crack_tip_x*1e3, self.crack_tip_y*1e3, 'ro', markersize=8, label='Crack Tip')
        
        # Refinement region
        circle = patches.Circle((self.crack_tip_x*1e3, self.crack_tip_y*1e3), 
                              self.tip_refinement_radius*1e3,
                              fill=False, edgecolor='green', linestyle='--', linewidth=2,
                              label='Refinement Zone')
        ax1.add_patch(circle)
        
        ax1.set_xlabel('x (mm)')
        ax1.set_ylabel('y (mm)')
        ax1.set_title('Robust Crack Mesh - Full Domain')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        
        # Plot 2: Zoomed view
        for triangle in triangles:
            triangle_nodes = nodes[triangle]
            triangle_closed = np.vstack([triangle_nodes, triangle_nodes[0]])
            ax2.plot(triangle_closed[:, 0]*1e3, triangle_closed[:, 1]*1e3, 'b-', linewidth=0.5, alpha=0.8)
        
        # Crack boundaries in zoom
        ax2.plot([0, self.crack_tip_x*1e3], [self.crack_y_bottom*1e3, self.crack_y_bottom*1e3], 'r-', linewidth=2)
        ax2.plot([0, self.crack_tip_x*1e3], [self.crack_y_top*1e3, self.crack_y_top*1e3], 'r-', linewidth=2)
        ax2.plot([self.crack_tip_x*1e3, self.crack_tip_x*1e3], [self.crack_y_bottom*1e3, self.crack_y_top*1e3], 'r-', linewidth=3)
        ax2.plot(self.crack_tip_x*1e3, self.crack_tip_y*1e3, 'ro', markersize=10, label='Crack Tip')
        
        # Set zoom window
        zoom_size = 0.15  # mm
        ax2.set_xlim((self.crack_tip_x - zoom_size*1e-3)*1e3, (self.crack_tip_x + zoom_size*1e-3)*1e3)
        ax2.set_ylim((self.crack_tip_y - zoom_size*1e-3)*1e3, (self.crack_tip_y + zoom_size*1e-3)*1e3)
        
        ax2.set_xlabel('x (mm)')
        ax2.set_ylabel('y (mm)')
        ax2.set_title('Crack Tip Refinement (Zoomed)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"\nMesh plot saved as {filename}")
        print(f"Mesh statistics:")
        print(f"  Nodes: {len(nodes)}")
        print(f"  Elements: {len(triangles)}")
        print(f"  Crack width: {self.crack_width*1e6:.1f} μm (consistent throughout)")
        print(f"  Crack tip location: ({self.crack_tip_x*1e3:.3f}, {self.crack_tip_y*1e3:.3f}) mm")

def main():
    """Main function"""
    generator = RobustCrackMesh()
    
    # Generate robust mesh
    nodes, triangles = generator.create_robust_mesh("robust_crack.msh")
    
    # Plot and analyze
    generator.plot_mesh(nodes, triangles, "robust_crack_mesh.png")
    
    print("\n✅ Robust mesh generation completed successfully!")
    print("✅ Mesh saved in MSH format for simulation import")

if __name__ == "__main__":
    main()