#pragma once
#include <vector>
#include <memory>
#include <cmath>
#include <Eigen/Dense>
#include "Node.h"

class Element {
public:
    virtual ~Element() = default;
    virtual Eigen::MatrixXd computeStiffnessMatrix() const = 0;
    virtual Eigen::VectorXd computeForceVector() const = 0;
    virtual std::vector<int> getGlobalDOFs() const = 0;
};

class Bar1D : public Element {
private:
    std::shared_ptr<Node> node1;
    std::shared_ptr<Node> node2;
    double E;
    double A;

public:
    Bar1D(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2, double E, double A)
        : node1(n1), node2(n2), E(E), A(A) {}

    double getLength() const {
        double dx = node2->x - node1->x;
        double dy = node2->y - node1->y;
        double dz = node2->z - node1->z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    Eigen::MatrixXd computeStiffnessMatrix() const override {
        double L = getLength();
        double k = (E * A) / L;
        Eigen::MatrixXd K(2, 2);
        K <<  k, -k,
             -k,  k;
        return K;
    }

    Eigen::VectorXd computeForceVector() const override {
        return Eigen::VectorXd::Zero(2);
    }

    std::vector<int> getGlobalDOFs() const override {
        return { node1->dofs[0], node2->dofs[0] };
    }
};