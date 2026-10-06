#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "Element.h"
#include "StressStrainPPE.h"

class Triangle3Node : public Element {
private:
    std::shared_ptr<Node> n1, n2, n3;
    Material material;
    double thickness;

public:
    Triangle3Node(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2, std::shared_ptr<Node> n3,
                  Material mat, double t = 1.0)
        : n1(n1), n2(n2), n3(n3), material(mat), thickness(t) {}

    // Strain-Displacement Matrix B (3x6 matrix for 2D 3-node triangle)
    Eigen::Matrix<double, 3, 6> computeBMatrix(double& area) const {
        double x1 = n1->x, y1 = n1->y;
        double x2 = n2->x, y2 = n2->y;
        double x3 = n3->x, y3 = n3->y;

        // Double area using determinant
        double detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        area = std::abs(detJ) / 2.0;

        double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
        double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

        Eigen::Matrix<double, 3, 6> B;
        B << b1,  0, b2,  0, b3,  0,
              0, c1,  0, c2,  0, c3,
             c1, b1, c2, b2, c3, b3;

        return B / (2.0 * area);
    }

    Eigen::MatrixXd computeStiffnessMatrix() const override {
        double area;
        auto B = computeBMatrix(area);
        auto D = material.getConstitutiveMatrix2D();

        // K_elem = thickness * area * B^T * D * B
        return thickness * area * B.transpose() * D * B;
    }

    Eigen::VectorXd computeForceVector() const override {
        return Eigen::VectorXd::Zero(6);
    }

    std::vector<int> getGlobalDOFs() const override {
        return { n1->dofs[0], n1->dofs[1],
                 n2->dofs[0], n2->dofs[1],
                 n3->dofs[0], n3->dofs[1] };
    }

    // Post-process element stresses given global displacement vector U
    StressResult computeElementStress(const Eigen::VectorXd& U_global) const {
        std::vector<int> dofs = getGlobalDOFs();
        Eigen::Vector<double, 6> u_elem;
        for (int i = 0; i < 6; ++i) {
            u_elem(i) = U_global(dofs[i]);
        }

        double area;
        auto B = computeBMatrix(area);
        Eigen::Vector3d strain = B * u_elem;

        return material.computeStress(strain);
    }

    std::vector<int> getNodeIDs() const override {
        return { n1->id, n2->id, n3->id };
    }
};
