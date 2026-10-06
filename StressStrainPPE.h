#pragma once
#include <Eigen/Dense>
#include <cmath>

struct StressResult {
    Eigen::Vector3d stressTensor; // [sigma_xx, sigma_yy, tau_xy]
    Eigen::Vector3d strainTensor; // [epsilon_xx, epsilon_yy, gamma_xy]
    double vonMises;
};

class Material {
public:
    double E;  // Young's Modulus
    double nu; // Poisson's Ratio

    Material(double E, double nu) : E(E), nu(nu) {}

    // Elastic Constitutive Matrix D (Plane Stress)
    Eigen::Matrix3d getConstitutiveMatrix2D() const {
        Eigen::Matrix3d D;
        double factor = E / (1.0 - nu * nu);
        D << 1.0,  nu, 0.0,
             nu,  1.0, 0.0,
             0.0, 0.0, (1.0 - nu) / 2.0;
        return factor * D;
    }

    // Post-processing helper: Compute Stress & Von Mises from Strain
    StressResult computeStress(const Eigen::Vector3d& strain) const {
        StressResult res;
        res.strainTensor = strain;
        res.stressTensor = getConstitutiveMatrix2D() * strain;

        double s_xx = res.stressTensor(0);
        double s_yy = res.stressTensor(1);
        double t_xy = res.stressTensor(2);

        // Von Mises stress for 2D Plane Stress
        res.vonMises = std::sqrt(s_xx * s_xx - s_xx * s_yy + s_yy * s_yy + 3.0 * t_xy * t_xy);
        return res;
    }
};
