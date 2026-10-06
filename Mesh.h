#pragma once
#include <vector>
#include <memory>
#include "Node.h"
#include "Element.h"

class Mesh {
public:
    std::vector<std::shared_ptr<Node>> nodes;
    std::vector<std::shared_ptr<Element>> elements;

    // Added default arguments = 0.0 for y and z
    std::shared_ptr<Node> addNode(double x, double y = 0.0, double z = 0.0) {
        int id = static_cast<int>(nodes.size());
        auto node = std::make_shared<Node>(id, x, y, z);
        node->dofs.push_back(id);
        nodes.push_back(node);
        return node;
    }

    void addElement(std::shared_ptr<Element> elem) {
        elements.push_back(elem);
    }

    int getTotalDOFs() const {
        return static_cast<int>(nodes.size());
    }
};