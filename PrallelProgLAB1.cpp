#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class GpuSpec {
public:
    std::string model;
    int memory; // in MB

    void Print() const {
        std::cout << "GPU Model: " << model << ", Memory: " << memory << " MB" << std::endl;
    }

    void Export(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::app);
        ofs << "GPU Model: " << model << ", Memory: " << memory << " MB" << std::endl;
        ofs.close();
    }

    void Import(const std::string& filename) {
        std::ifstream ifs(filename);
        if (ifs.is_open()) {
            std::getline(ifs, model, ',');
            ifs >> memory;
            ifs.close();
        }
    }
};

class CpuSpec {
public:
    std::string model;
    int cores;

    void Print() const {
        std::cout << "CPU Model: " << model << ", Cores: " << cores << std::endl;
    }

    void Export(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::app);
        ofs << "CPU Model: " << model << ", Cores: " << cores << std::endl;
        ofs.close();
    }

    void Import(const std::string& filename) {
        std::ifstream ifs(filename);
        if (ifs.is_open()) {
            std::getline(ifs, model, ',');
            ifs >> cores;
            ifs.close();
        }
    }
};

class RamSpec {
public:
    int size; // in GB

    void Print() const {
        std::cout << "RAM Size: " << size << " GB" << std::endl;
    }

    void Export(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::app);
        ofs << "RAM Size: " << size << " GB" << std::endl;
        ofs.close();
    }

    void Import(const std::string& filename) {
        std::ifstream ifs(filename);
        if (ifs.is_open()) {
            ifs >> size;
            ifs.close();
        }
    }
};

class LanSpec {
public:
    std::string model;

    void Print() const {
        std::cout << "LAN Model: " << model << std::endl;
    }

    void Export(const std::string& filename) const {
        std::ofstream ofs(filename, std::ios::app);
        ofs << "LAN Model: " << model << std::endl;
        ofs.close();
    }

    void Import(const std::string& filename) {
        std::ifstream ifs(filename);
        if (ifs.is_open()) {
            std::getline(ifs, model);
            ifs.close();
        }
    }
};

class ClusterNode {
public:
    GpuSpec gpu;
    CpuSpec cpu;
    RamSpec ram;
    LanSpec lan;

    void Print() const {
        gpu.Print();
        cpu.Print();
        ram.Print();
        lan.Print();
    }

    void Export(const std::string& filename) const {
        gpu.Export(filename);
        cpu.Export(filename);
        ram.Export(filename);
        lan.Export(filename);
    }

    void Import(const std::string& filename) {
        gpu.Import(filename);
        cpu.Import(filename);
        ram.Import(filename);
        lan.Import(filename);
    }
};

class Cluster {
public:
    std::vector<ClusterNode> nodes;

    void AddNode(const ClusterNode& node) {
        nodes.push_back(node);
    }

    void Print() const {
        for (const auto& node : nodes) {
            node.Print();
            std::cout << "------------------------" << std::endl;
        }
    }

    void Export(const std::string& filename) const {
        for (const auto& node : nodes) {
            node.Export(filename);
            // Add a separator for clarity
            std::ofstream ofs(filename, std::ios::app);
            ofs << "------------------------" << std::endl;
            ofs.close();
        }
    }

    void Import(const std::string& filename) {
        // This is a simplified import. In a real scenario, you would need to handle multiple nodes.
        ClusterNode node;
        node.Import(filename); // Import the last node
        AddNode(node);
    }
};

#include <gtest/gtest.h>

TEST(GpuSpecTest, PrintAndExport) {
    GpuSpec gpu;
    gpu.model = "NVIDIA RTX 3080";
    gpu.memory = 10240;

    testing::internal::CaptureStdout();
    gpu.Print();
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("NVIDIA RTX 3080"), std::string::npos);

    gpu.Export("gpu_test.txt");
}

TEST(CpuSpecTest, PrintAndExport) {
    CpuSpec cpu;
    cpu.model = "Intel i9-10900K";
    cpu.cores = 10;

    testing::internal::CaptureStdout();
    cpu.Print();
    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("Intel i9-10900K"), std::string::npos);

    cpu.Export("cpu_test.txt");
}

// Добавьте аналогичные тесты для RamSpec и LanSpec.

TEST(ClusterNodeTest, PrintAndExport) {
    ClusterNode node;

    node.gpu.model = "NVIDIA RTX 3080";
    node.gpu.memory = 10240;

    node.cpu.model = "Intel i9-10900K";
    node.cpu.cores = 10;

    node.ram.size = 32; // in GB

    node.lan.model = "Realtek RTL8111";

    testing::internal::CaptureStdout();
    node.Print();

    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("NVIDIA RTX 3080"), std::string::npos);

    node.Export("node_test.txt");
}

TEST(ClusterTest, AddNodeAndPrint) {
    Cluster cluster;

    ClusterNode node1;

    node1.gpu.model = "NVIDIA RTX 3080";
    node1.cpu.model = "Intel i9-10900K";

    cluster.AddNode(node1);

    testing::internal::CaptureStdout();
    cluster.Print();

    std::string output = testing::internal::GetCapturedStdout();

    EXPECT_NE(output.find("NVIDIA RTX 3080"), std::string::npos);
}
