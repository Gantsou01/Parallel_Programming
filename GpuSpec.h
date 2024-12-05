#ifndef GPUSPEC_H
#define GPUSPEC_H

#include <string>
#include <iostream>
#include <fstream>

class GpuSpec {
private:
    std::string model;
    int memory;  // Б ла

public:
    GpuSpec(const std::string& model = "", int memory = 0)
        : model(model), memory(memory) {}

    void Print() const {
        std::cout << "GPU Model: " << model << ", Memory: " << memory << " MB\n";
    }

    void Import(const std::string& filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            file >> model >> memory;
        }
    }

    void Export(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << model << " " << memory << "\n";
        }
    }
};

#endif

