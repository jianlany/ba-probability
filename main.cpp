#include <iostream>
#include <vector>
#include <string>
#include <fstream>

struct BA {
    double d1, d2;
    double q;
};
using BA_data = std::vector<BA>;

BA_data read_BA_file(std::string path);

int main() {

    std::string ba_path = "/home/joswald1/development/ba-probability"
                          "/ba-PE_50chains_160beads-1_EEE.txt";
   
    auto ba_data = read_BA_file(ba_path); 
}


BA_data read_BA_file(std::string path) {
    BA_data ba_data;
    std::fstream fid(path);

    while (fid) {
        BA ba;
        fid >> ba.d1 >> ba.d2 >> ba.q;
        if (fid) {
            ba_data.push_back(ba);
        }
    }

    return ba_data;
}

