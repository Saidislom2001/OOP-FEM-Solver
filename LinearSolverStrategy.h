#pragma once
#include <Eigen/Sparse>
#include <memory>

class LinearSolverStrategy {
public:
    virtual ~LinearSolverStrategy() = default;
    virtual Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& F) = 0;
};

// Direct Sparse Cholesky Solver (For Symmetric Positive Definite systems)
class CholeskySolver : public LinearSolverStrategy {
public:
    Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& F) override {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(K);
        return solver.solve(F);
    }
};

// Conjugate Gradient Solver (For large-scale systems)
class CGSolver : public LinearSolverStrategy {
public:
    Eigen::VectorXd solve(const Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& F) override {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
        cg.compute(K);
        return cg.solve(F);
    }
};
