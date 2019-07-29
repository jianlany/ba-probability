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

    // Bounds on bond-length histogram - TODO: should be set by user.
    const int M=100, N=100;
    const auto xlo = 2.0,  xhi = 2.9;
    const auto ylo = 70.0, yhi = 180.0;
    // Compute grid spacing.
    const auto dx = (xhi-xlo) / (M-1);
    const auto dy = (yhi-ylo) / (N-1);
    // kernel width is 1.5 grid spacings.
    const auto wx = 1.5*dx;
    const auto wy = 1.5*dy;
    const auto Z = 1.0 / (2.0 * acos(-1.0) * sqrt(wx*wy)) / (2.0 * ba_data.size());
    double phi[M*N] = {0.0};
    #pragma acc data copy(data[:3*n_samples], phi)
    {
    #pragma acc parallel loop collapse(2) 
    for (int j=0; j<M; ++j) {
        for (int k=0; k<N; ++k) {
            for (int i=0; i<n_samples; ++i) {
                auto d1 = data[3*i];
                auto d2 = data[3*i+1];
                auto q = data[3*i+2];
                auto xj = xlo + dx*j;
                auto yk = ylo + dy*k;
                auto x1 = (d1-xj) / wx;
                auto x2 = (d2-xj) / wx;
                auto y = (q-yk) / wy;
                phi[k + j*N] += Z*exp(-0.5*(x1*x1 + y*y));
                phi[k + j*N] += Z*exp(-0.5*(x2*x2 + y*y));
            }
        }
    }
    }
    auto sum = 0.0;
    for (int i=0; i<M-1; ++i) {
        for (int j=0; j<N-1; ++j) {
            auto p1 = phi[j   + i*N];
            auto p2 = phi[j   + (i+1)*N];
            auto p3 = phi[j+1 + (i+1)*N];
            auto p4 = phi[j+1 + i*N];
            sum += 0.25*(p1+p2+p3+p4) * dx * dy;
        }
    }
    // Sum should be close to 1.0.
    std::cout << "Sum = " << sum << "\n";
    std::fstream fid("p.txt", std::ios::out);
    for (int i=0; i<M; ++i) {
        if (i > 0) fid << "\n";
        for (int j=0; j<N; ++j) {
            if (j > 0) fid << " ";
            fid << phi[j + i*N];
        }
    }
}


// Reads the BA datafile output from CG-distributions.  Each row of the file 
// should contain: [bond_len1, bond_len2, bond_angle].
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

