#pragma once
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "Mesh.h"

class Plotter {
public:
    static void plot1D(const Mesh& mesh, const Eigen::VectorXd& U, const std::string& outputImage = "displacement.png") {
        // Export CSV data
        std::ofstream dataFile("results.csv");
        dataFile << "x,u\n";
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            dataFile << mesh.nodes[i]->x << "," << U(i) << "\n";
        }
        dataFile.close();

        // Write lightweight Python script
        std::ofstream pyFile("plot_results.py");
        pyFile << "import pandas as pd\n";
        pyFile << "import matplotlib.pyplot as plt\n\n";
        pyFile << "df = pd.read_csv('results.csv')\n";
        pyFile << "plt.figure(figsize=(8, 5))\n";
        pyFile << "plt.plot(df['x'], df['u'], 'o-', label='Displacement (m)', color='b', linewidth=2)\n";
        pyFile << "plt.title('FEM 1D Bar Displacement')\n";
        pyFile << "plt.xlabel('Position X (m)')\n";
        pyFile << "plt.ylabel('Displacement U (m)')\n";
        pyFile << "plt.grid(True)\n";
        pyFile << "plt.legend()\n";
        pyFile << "plt.savefig('" << outputImage << "')\n";
        pyFile << "print('Plot saved as " << outputImage << "')\n";
        pyFile.close();

        // Trigger Python script execution automatically
        int ret = std::system("python3 plot_results.py");
        if (ret != 0) {
            std::cout << "Data saved to 'results.csv'. Run 'python3 plot_results.py' manually to display the plot.\n";
        }
    }
};
