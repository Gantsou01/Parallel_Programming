#ifndef CPUSPEC_H
#define CPUSPEC_H

#include <string>
#include <iostream>
#include <fstream>

class CpuSpec {
private:
    std::string model;
    int cores;

public:
    CpuSpec(const std::string& model = "", int cores = 0)
        : model(model), cores(cores) {}

    void Print() const {
        std::cout << "CPU Model: " << model << ", Cores: " << cores << "\n";
    }

    void Import(const std::string& filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            file >> model >> cores;
        }
    }

    void Export(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << model << " " << cores << "\n";
        }
    }
};

#endif
