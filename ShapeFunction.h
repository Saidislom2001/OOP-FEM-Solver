#pragma once
#include <Eigen/Dense>

// Abstract Shape Function class for Isoparametric Elements
class ShapeFunction2D {
public:
    virtual ~ShapeFunction2D() = default;
    
    // Evaluate N_i at parametric coordinates (xi, eta)
    virtual Eigen::VectorXd evaluateN(double xi, double eta) const = 0;
    
    // Evaluate derivatives [dN/dxi, dN/deta]
    virtual Eigen::MatrixXd evaluateDerivatives(double xi, double eta) const = 0;
};

// 4-Node Bilinear Quadrilateral Element (Quad4)
class Quad4Shape : public ShapeFunction2D {
public:
    Eigen::VectorXd evaluateN(double xi, double eta) const override {
        Eigen::Vector4d N;
        N(0) = 0.25 * (1 - xi) * (1 - eta);
        N(1) = 0.25 * (1 + xi) * (1 - eta);
        N(2) = 0.25 * (1 + xi) * (1 + eta);
        N(3) = 0.25 * (1 - xi) * (1 + eta);
        return N;
    }

    Eigen::MatrixXd evaluateDerivatives(double xi, double eta) const override {
        Eigen::MatrixXd dN(2, 4);
        // dN/dxi
        dN(0, 0) = -0.25 * (1 - eta);
        dN(0, 1) =  0.25 * (1 - eta);
        dN(0, 2) =  0.25 * (1 + eta);
        dN(0, 3) = -0.25 * (1 + eta);
        // dN/deta
        dN(1, 0) = -0.25 * (1 - xi);
        dN(1, 1) = -0.25 * (1 + xi);
        dN(1, 2) =  0.25 * (1 + xi);
        dN(1, 3) =  0.25 * (1 - xi);
        return dN;
    }
};
