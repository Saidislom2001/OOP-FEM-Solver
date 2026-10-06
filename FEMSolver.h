#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "Mesh.h"

struct DirichletBC {
    int nodeID;
    int dof;
    double value;
};

struct NeumannBC {
    int nodeID;
    int dof;
    double force;
};

class FEMSolver {
private:
    Mesh& mesh;
    Eigen::MatrixXd K_global;
    Eigen::VectorXd F_global;
    Eigen::VectorXd U;

    std::vector<DirichletBC> dirichletBCs;
    std::vector<NeumannBC> neumannBCs;

public:
    explicit FEMSolver(Mesh& mesh) : mesh(mesh) {}

    void addDirichletBC(int nodeID, int dof, double value) {
        dirichletBCs.push_back({nodeID, dof, value});
    }

    void addNeumannBC(int nodeID, int dof, double force) {
        neumannBCs.push_back({nodeID, dof, force});
    }

    void assemble() {
        int nDOFs = mesh.getTotalDOFs();
        K_global = Eigen::MatrixXd::Zero(nDOFs, nDOFs);
        F_global = Eigen::VectorXd::Zero(nDOFs);

        for (const auto& elem : mesh.elements) {
            Eigen::MatrixXd K_elem = elem->computeStiffnessMatrix();
            Eigen::VectorXd F_elem = elem->computeForceVector();
            std::vector<int> dofs = elem->getGlobalDOFs();

            for (size_t i = 0; i < dofs.size(); ++i) {
                F_global(dofs[i]) += F_elem(i);
                for (size_t j = 0; j < dofs.size(); ++j) {
                    K_global(dofs[i], dofs[j]) += K_elem(i, j);
                }
            }
        }

        for (const auto& nbc : neumannBCs) {
            int globalDOF = mesh.nodes[nbc.nodeID]->dofs[nbc.dof];
            F_global(globalDOF) += nbc.force;
        }
    }

    void applyDirichletBCs() {
        for (const auto& dbc : dirichletBCs) {
            int dof = mesh.nodes[dbc.nodeID]->dofs[dbc.dof];

            K_global.row(dof).setZero();
            K_global.col(dof).setZero();
            
            K_global(dof, dof) = 1.0;
            F_global(dof) = dbc.value;
        }
    }

    void solve() {
        assemble();
        applyDirichletBCs();
        U = K_global.colPivHouseholderQr().solve(F_global);
    }

    const Eigen::VectorXd& getDisplacements() const { 
    return U; 
    }

    void printResults() const {
        std::cout << "\n--- Displacements (U) ---\n";
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            std::cout << "Node " << i << ": " << U(i) << " m\n";
        }
    }
};