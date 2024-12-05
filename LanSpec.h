#ifndef LANSPEC_H
#define LANSPEC_H

#include <string>
#include <iostream>
#include <fstream>

class LanSpec {
private:
    std::string type;  // например, Ethernet
    int speed;         // скорость в ћбит/с

public:
    LanSpec(const std::string& type = "", int speed = 0)
        : type(type), speed(speed) {}

    void Print() const {
        std::cout << "LAN Type: " << type << ", Speed: " << speed << " Mbps\n";
    }

    void Import(const std::string& filename) {
        std::ifstream file(filename);
        if (file.is_open()) {
            file >> type >> speed;
        }
    }

    void Export(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << type << " " << speed << "\n";
        }
    }
};

#endif
