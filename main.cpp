#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

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
    const double *data = &ba_data[0].d1;
    int n_samples = ba_data.size();

    const int M=100, N=100;
    double phi[M*N] = {0.0};
    #pragma acc data copy(data[:3*n_samples], phi)
    {
    #pragma acc parallel loop collapse(2) 
    for (int j=0; j<M; ++j) {
        for (int k=0; k<N; ++k) {
            for (int i=0; i<n_samples; ++i) {
                double d = data[3*i];
                auto x = (d-2.2) / 0.05;
                phi[k + j*N] += exp(-x);
            }
        }
    }
    }
    auto sum = 0.0;
    for (int i=0; i<M*N; ++i) sum += phi[i];
    std::cout << "sum = " << sum << "\n";
}


BA_data read_BA_file(std::string path) {
    BA_data ba_data;
    std::ifstream fid(path, std::ios::in);
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

