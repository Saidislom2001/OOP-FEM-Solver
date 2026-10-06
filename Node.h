#pragma once
#include <vector>

class Node {
public:
    int id;
    double x, y, z;
    std::vector<int> dofs;

    Node(int id, double x, double y = 0.0, double z = 0.0) 
        : id(id), x(x), y(y), z(z) {}
};

    