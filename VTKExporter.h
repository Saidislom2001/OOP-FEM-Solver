#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "Mesh.h"

class VTKExporter {
public:
    static void exportToVTK(const std::string& filename, const Mesh& mesh, const Eigen::VectorXd& U) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing.\n";
            return;
        }

        // 1. VTK Header
        file << "# vtk DataFile Version 3.0\n";
        file << "FEM Solver Results\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n\n";

        // 2. Node Coordinates
        file << "POINTS " << mesh.nodes.size() << " double\n";
        for (const auto& node : mesh.nodes) {
            file << node->x << " " << node->y << " " << node->z << "\n";
        }
        file << "\n";

        // 3. Elements (Cells)
        // For 1D bar elements (VTK_LINE = 3), each cell has 2 points
        file << "CELLS " << mesh.elements.size() << " " << mesh.elements.size() * 3 << "\n";
        for (const auto& elem : mesh.elements) {
            auto dofs = elem->getGlobalDOFs();
            file << "2 " << dofs[0] << " " << dofs[1] << "\n";
        }
        file << "\n";

        // Cell Types (3 = VTK_LINE)
        file << "CELL_TYPES " << mesh.elements.size() << "\n";
        for (size_t i = 0; i < mesh.elements.size(); ++i) {
            file << "3\n";
        }
        file << "\n";

        // 4. Nodal Solution Data (Displacements)
        file << "POINT_DATA " << mesh.nodes.size() << "\n";
        file << "SCALARS Displacement double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            file << U(i) << "\n";
        }

        file.close();
        std::cout << "VTK file written successfully to: " << filename << "\n";
    }
};
