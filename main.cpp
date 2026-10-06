#include <iostream>
#include <memory>
#include "Mesh.h"
#include "FEMSolver.h"
#include "Triangle3Node.h"
#include "StressStrainPPE.h"
#include "LinearSolverStrategy.h" // <--- Include added here
#include "VTKExporter.h"

int main() {
    // 1. Define Mesh Parameters
    Mesh mesh;
    double Length = 2.0; // Beam length in meters
    double Height = 0.5; // Beam height in meters
    int nx = 4;          // Elements along length
    int ny = 2;          // Elements along height

    double dx = Length / nx;
    double dy = Height / ny;

    // Create 2D Nodes
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            mesh.addNode2D(i * dx, j * dy);
        }
    }

    auto getNodeIndex = [nx](int i, int j) {
        return j * (nx + 1) + i;
    };

    // 2. Define Material Properties
    Material steel(210e9, 0.3); // E = 210 GPa, nu = 0.3
    double thickness = 0.01;    // 10 mm plate thickness

    // 3. Create 2D CST Elements
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            auto n0 = mesh.nodes[getNodeIndex(i, j)];
            auto n1 = mesh.nodes[getNodeIndex(i + 1, j)];
            auto n2 = mesh.nodes[getNodeIndex(i + 1, j + 1)];
            auto n3 = mesh.nodes[getNodeIndex(i, j + 1)];

            mesh.addElement(std::make_shared<Triangle3Node>(n0, n1, n2, steel, thickness));
            mesh.addElement(std::make_shared<Triangle3Node>(n0, n2, n3, steel, thickness));
        }
    }

    // 4. Set Boundary Conditions
    FEMSolver solver(mesh, std::make_unique<CholeskySolver>());

    // Fix Left Edge (X = 0)
    for (int j = 0; j <= ny; ++j) {
        int fixedNodeID = getNodeIndex(0, j);
        solver.addDirichletBC(fixedNodeID, 0, 0.0); // Ux = 0
        solver.addDirichletBC(fixedNodeID, 1, 0.0); // Uy = 0
    }

    // Apply Downward Load at Top-Right Corner Node
    int loadedNodeID = getNodeIndex(nx, ny);
    double downwardForce = -10000.0; // -10 kN in Y direction
    solver.addNeumannBC(loadedNodeID, 1, downwardForce);

    // 5. Assemble and Solve
    std::cout << "Assembling and solving 2D FEM System (" << mesh.getTotalDOFs() << " DOFs)...\n";
    solver.solve();

    // 6. Post-Processing: Compute Stresses
    const auto& U = solver.getDisplacements();
    std::cout << "\n--- Element Stresses (Von Mises) ---\n";
    
    for (size_t i = 0; i < mesh.elements.size(); ++i) {
        auto triElem = std::dynamic_pointer_cast<Triangle3Node>(mesh.elements[i]);
        if (triElem) {
            StressResult stress = triElem->computeElementStress(U);
            std::cout << "Element " << i 
                      << " | Von Mises Stress: " << stress.vonMises / 1e6 << " MPa"
                      << " | Sigma_XX: " << stress.stressTensor(0) / 1e6 << " MPa\n";
        }
    }

    // 7. Export Results to VTK
    VTKExporter::exportToVTK("cantilever_2d.vtk", mesh, U);

    return 0;
}