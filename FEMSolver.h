#pragma once
#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Mesh.h"
#include "LinearSolverStrategy.h"

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
    Eigen::SparseMatrix<double> K_sparse;
    Eigen::VectorXd F_global;
    Eigen::VectorXd U;

    std::vector<DirichletBC> dirichletBCs;
    std::vector<NeumannBC> neumannBCs;
    std::unique_ptr<LinearSolverStrategy> strategy;

public:
    // Constructor accepts an optional sparse solver strategy
    explicit FEMSolver(Mesh& mesh, std::unique_ptr<LinearSolverStrategy> solverStrategy = nullptr) 
        : mesh(mesh), strategy(std::move(solverStrategy)) {}

    void addDirichletBC(int nodeID, int dof, double value) {
        dirichletBCs.push_back({nodeID, dof, value});
    }

    void addNeumannBC(int nodeID, int dof, double force) {
        neumannBCs.push_back({nodeID, dof, force});
    }

    // --- DENSE SOLVER METHODS (Restored) ---
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

    // --- SPARSE SOLVER METHODS ---
    void assembleSparse() {
        int nDOFs = mesh.getTotalDOFs();
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(mesh.elements.size() * 36); 
        F_global = Eigen::VectorXd::Zero(nDOFs);

        for (const auto& elem : mesh.elements) {
            Eigen::MatrixXd K_elem = elem->computeStiffnessMatrix();
            Eigen::VectorXd F_elem = elem->computeForceVector();
            std::vector<int> dofs = elem->getGlobalDOFs();

            for (size_t i = 0; i < dofs.size(); ++i) {
                F_global(dofs[i]) += F_elem(i);
                for (size_t j = 0; j < dofs.size(); ++j) {
                    tripletList.push_back(Eigen::Triplet<double>(dofs[i], dofs[j], K_elem(i, j)));
                }
            }
        }

        for (const auto& nbc : neumannBCs) {
            int globalDOF = mesh.nodes[nbc.nodeID]->dofs[nbc.dof];
            F_global(globalDOF) += nbc.force;
        }

        double penalty = 1e12; 
        for (const auto& dbc : dirichletBCs) {
            int dof = mesh.nodes[dbc.nodeID]->dofs[dbc.dof];
            tripletList.push_back(Eigen::Triplet<double>(dof, dof, penalty));
            F_global(dof) += dbc.value * penalty;
        }

        K_sparse.resize(nDOFs, nDOFs);
        K_sparse.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void solve() {
        if (strategy) {
            assembleSparse();
            U = strategy->solve(K_sparse, F_global);
        } else {
            assemble();
            applyDirichletBCs();
            U = K_global.colPivHouseholderQr().solve(F_global);
        }
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