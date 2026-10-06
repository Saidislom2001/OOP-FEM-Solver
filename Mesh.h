#pragma once
#include <vector>
#include <memory>
#include "Node.h"
#include "Element.h"

class Mesh {
public:
    std::vector<std::shared_ptr<Node>> nodes;
    std::vector<std::shared_ptr<Element>> elements;

    // 1D Node creation (1 DOF)
    std::shared_ptr<Node> addNode(double x, double y = 0.0, double z = 0.0) {
        int id = static_cast<int>(nodes.size());
        auto node = std::make_shared<Node>(id, x, y, z);
        node->dofs.push_back(id);
        nodes.push_back(node);
        return node;
    }

    // 2D Node creation (2 DOFs: Ux, Uy)
    std::shared_ptr<Node> addNode2D(double x, double y) {
        int id = static_cast<int>(nodes.size());
        auto node = std::make_shared<Node>(id, x, y, 0.0);
        int startDOF = id * 2;
        node->dofs.push_back(startDOF);     // Ux
        node->dofs.push_back(startDOF + 1); // Uy
        nodes.push_back(node);
        return node;
    }

    void addElement(std::shared_ptr<Element> elem) {
        elements.push_back(elem);
    }

    // Single getTotalDOFs implementation
    int getTotalDOFs() const {
        int total = 0;
        for (const auto& node : nodes) {
            total += static_cast<int>(node->dofs.size());
        }
        return total;
    }
};