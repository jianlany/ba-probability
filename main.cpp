#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

struct BA {
    double d1, d2;
    double q;
};
using BA_data = std::vector<BA>;

BA_data read_BA_file(std::string path);


int main() {
    std::string ba_path = "/home/joswald1/development/ba-probability"
                          "/ba-PE_50chains_160beads-1_EEE.txt";
    double sum = 0.0;
    auto ba_data = read_BA_file(ba_path); 
    for (int i=0; i<ba_data.size(); ++i) {
        #pragma acc parallel loop collapse(2)
        for (int j=0; j<20; ++j) {
            for (int k=0; k<20; ++k) {
               sum += exp(-(ba_data[i].d1-2.2) / 0.05);
            }
        }
    }
    std::cout << "sum = " << sum << "\n";
}


BA_data read_BA_file(std::string path) {
    BA_data ba_data;
    std::fstream fid(path);
    if (!fid) {
        std::cerr << "Could not open " << path << ".\n";
        exit(1);
    }
    while (fid) {
        BA ba;
        fid >> ba.d1 >> ba.d2 >> ba.q;
        if (fid) {
            ba_data.push_back(ba);
        }
    }
    return ba_data;
}

