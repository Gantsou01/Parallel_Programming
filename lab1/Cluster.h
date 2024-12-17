#ifndef CLUSTER_H
#define CLUSTER_H

#include "ClusterNode.h"
#include <vector>
#include <iostream>

class Cluster {
private:
    std::vector<ClusterNode> nodes;

public:
    void AddNode(const ClusterNode& node) {
        nodes.push_back(node);
    }

    void Print() const {
        std::cout << "Cluster Overview:\n";
        for (const auto& node : nodes) {
            node.Print();
        }
    }

    void Import(const std::string& filename) {
        for (auto& node : nodes) {
            node.Import(filename);
        }
    }

    void Export(const std::string& filename) const {
        for (const auto& node : nodes) {
            node.Export(filename);
        }
    }
};

#endif

