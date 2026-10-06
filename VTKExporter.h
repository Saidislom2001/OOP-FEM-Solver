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
            std::cerr << "Error: Could not open file " << filename << "\n";
            return;
        }

        file << "# vtk DataFile Version 3.0\nFEM Solver Results\nASCII\nDATASET UNSTRUCTURED_GRID\n\n";

        // Nodes
        file << "POINTS " << mesh.nodes.size() << " double\n";
        for (const auto& node : mesh.nodes) {
            file << node->x << " " << node->y << " " << node->z << "\n";
        }
        file << "\n";

        // Dynamic Elements
        int totalCellNodes = 0;
        for (const auto& elem : mesh.elements) { 
            totalCellNodes += elem->getNodeIDs().size(); 
        }
        
        file << "CELLS " << mesh.elements.size() << " " << mesh.elements.size() + totalCellNodes << "\n";
        for (const auto& elem : mesh.elements) {
            auto nodes = elem->getNodeIDs();
            file << nodes.size();
            for (int n : nodes) file << " " << n;
            file << "\n";
        }
        file << "\n";

        // Dynamic Cell Types
        file << "CELL_TYPES " << mesh.elements.size() << "\n";
        for (const auto& elem : mesh.elements) {
            size_t nNodes = elem->getNodeIDs().size();
            if (nNodes == 2) file << "3\n";      // VTK_LINE
            else if (nNodes == 3) file << "5\n"; // VTK_TRIANGLE
            else if (nNodes == 4) file << "9\n"; // VTK_QUAD
            else file << "1\n";                  // VTK_VERTEX
        }
        file << "\n";

        // Displacements (Exported as VECTORS to support 2D/3D ParaView Glyphs)
        file << "POINT_DATA " << mesh.nodes.size() << "\n";
        file << "VECTORS Displacement double\n";
        for (const auto& node : mesh.nodes) {
            double u = (node->dofs.size() > 0) ? U(node->dofs[0]) : 0.0;
            double v = (node->dofs.size() > 1) ? U(node->dofs[1]) : 0.0;
            double w = (node->dofs.size() > 2) ? U(node->dofs[2]) : 0.0;
            file << u << " " << v << " " << w << "\n";
        }

        file.close();
        std::cout << "VTK file written successfully to: " << filename << "\n";
    }
};