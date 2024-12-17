#ifndef CLUSTERNODE_H
#define CLUSTERNODE_H

#include "GpuSpec.h"
#include "CpuSpec.h"
#include "RpmSpec.h"
#include "LanSpec.h"
#include <vector>
#include <iostream>

class ClusterNode {
private:
    GpuSpec gpu;
    CpuSpec cpu;
    RpmSpec rpm;
    LanSpec lan;

public:
    ClusterNode(const GpuSpec& gpu, const CpuSpec& cpu, const RpmSpec& rpm, const LanSpec& lan)
        : gpu(gpu), cpu(cpu), rpm(rpm), lan(lan) {}

    void Print() const {
        std::cout << "Cluster Node Specifications:\n";
        gpu.Print();
        cpu.Print();
        rpm.Print();
        lan.Print();
    }

    void Import(const std::string& filename) {
        gpu.Import(filename);
        cpu.Import(filename);
        rpm.Import(filename);
        lan.Import(filename);
    }

    void Export(const std::string& filename) const {
        gpu.Export(filename);
        cpu.Export(filename);
        rpm.Export(filename);
        lan.Export(filename);
    }
};

#endif

