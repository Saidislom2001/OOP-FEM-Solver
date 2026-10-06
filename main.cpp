#include "Mesh.h"
#include "FEMSolver.h"
#include "VTKExporter.h"
#include "Plotter.h"

int main() {
    Mesh mesh;
    // Discretize a 2.0m bar into 10 elements
    int numElements = 10;
    double L_total = 2.0;
    double dx = L_total / numElements;

    for (int i = 0; i <= numElements; ++i) {
        mesh.addNode(i * dx);
    }

    double E = 210e9; // 210 GPa
    double A = 0.01;  // 0.01 m^2

    for (int i = 0; i < numElements; ++i) {
        mesh.addElement(std::make_shared<Bar1D>(mesh.nodes[i], mesh.nodes[i+1], E, A));
    }

    FEMSolver solver(mesh);
    solver.addDirichletBC(0, 0, 0.0);           // Fixed at X=0
    solver.addNeumannBC(numElements, 0, 10000.0); // 10 kN force at end

    solver.solve();
    solver.printResults();

    // 1. Export results for ParaView (VTK)
    VTKExporter::exportToVTK("results.vtk", mesh, solver.getDisplacements());

    // 2. Generate 1D plot image via Python
    Plotter::plot1D(mesh, solver.getDisplacements(), "displacement_plot.png");

    return 0;
}