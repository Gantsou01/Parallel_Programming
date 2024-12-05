#ifndef RPMSPEC_H
#define RPMSPEC_H

#include <string>
#include <iostream>
#include <fstream>

class RpmSpec {
private:
    std::string type;  // HDD ่๋่ SSD
    int capacity;      // โ รม

public:
    RpmSpec(const std::string& type = "", int capacity = 0)
        : type(type), capacity(capacity) {}

    void Print() const {
        std::cout << "RPM Type: " << type << ", Capacity: " << capacity << " GB\n";
    }

    void Import(const std::string& filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            file >> type >> capacity;
        }
    }

    void Export(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << type << " " << capacity << "\n";
        }
    }
};

#endif

