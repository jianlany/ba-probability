#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <fstream>
#include <math.h>
const double pi = 3.14159265358979323846;
const std::string help_msg = 
"Usage:\n"
    "\tba-probability [path] (options).\n"
    "\toptions:\n"
        "\t\t--num_theta: number of bins on angle direction. Default: 100\n"
        "\t\t--num_length: number of bins on bond direction. Default: 100\n"
        "\t\t--d_width_factor: number of bins that the gaussion distribution will spread over on bond direction. Default: 1.5\n"
        "\t\t--theta_width_factor: same but on angle direction. Default: 1.5\n"
        "\t\t--d_range: the upper and lower limit of bond. Default: 2.0  2.9\n"
        "\t\t--theta_range: the upper and lower limit of angle. Default: 70  180\n"
        "\t\t--help: print this message.\n";
struct BA {
    double d1, d2;
    double q;
};
using BA_data = std::vector<BA>;

BA_data read_BA_file(std::string path);


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage:\n  ./ba-probability [path] (options)\n.";
        return 1;
    }
    std::string ba_path = argv[1];

    // Bounds on bond-length histogram - TODO: should be set by user.
    int M=100, N=100;
    auto xlo = 2.0,  xhi = 2.9;
    auto ylo = 70.0, yhi = 180.0;
    // How many grid spacings the distribution is spread over.
    double d_width_factor = 1.5;
    double q_width_factor = 1.5;
    // Loop over remaining arguments and set optional parameters.
    for (int i=2; i<argc; ++i) {
        if (!strcmp(argv[i], "--help")) {
            std::cout << help_msg;
            exit(0);
        }
        else if (!strcmp(argv[i], "--num_theta") && ++i < argc) {
            N = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--num_length") && ++i < argc) {
            M = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--d_width_factor") && ++i < argc) {
            d_width_factor = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--q_width_factor") && ++i < argc) {
            q_width_factor = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--d_range") && (i+2) < argc) {
            xlo = atoi(argv[++i]);
            xhi = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--theta_range") && (i+2) < argc) {
            ylo = atoi(argv[++i]);
            yhi = atoi(argv[++i]);
        }
        else{
            std::cerr << "Invalid input arguments: " << argv[i] << " \n";
            std::cerr << help_msg;
            exit(1);
        }

    }

    // Compute grid spacing.
    const auto dx = (xhi-xlo) / (M-1);
    const auto dy = (yhi-ylo) / (N-1);
    // kernel width is 1.5 grid spacings.
    const auto wx = d_width_factor*dx;
    const auto wy = q_width_factor*dy;
    auto ba_data = read_BA_file(ba_path);
    const double *data = &ba_data[0].d1;
    int n_samples = ba_data.size();
    const auto Z = 1.0 / (2.0 * pi * wx*wy) / (2.0 * ba_data.size());
    double *phi = new double[M*N];
    double *phi_l = new double[M*N];
    double *phi_q = new double[M*N];
    double *phi_lq = new double[M*N];
    for (int i=0; i<M*N; ++i) {
        phi[i] = 0.0;
        phi_l[i] = 0.0;
        phi_q[i] = 0.0;
        phi_lq[i] = 0.0;
    }
    #pragma acc data copy(data[:3*n_samples], phi[:M*N], phi_l[:M*N], phi_q[:M*N], phi_lq[:M*N])
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
                auto x1 = (xj-d1) / wx;
                auto x2 = (xj-d2) / wx;
                auto y = (yk-q) / wy;
                auto exp1 = exp(-0.5*(x1*x1 + y*y));
                auto exp2 = exp(-0.5*(x2*x2 + y*y));
                phi[k + j*N] += Z*exp1 + Z*exp2;
                phi_l[k + j*N] += Z*exp1*(-x1/wx) + Z*exp2*(-x2/wx);
                phi_q[k + j*N] += Z*exp1*(-y/wy) + Z*exp2*(-y/wy);
                phi_lq[k + j*N] += Z*exp1*(y/wy)*(x1/wx) +
                                   Z*exp2*(y/wy)*(x2/wx);
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
    std::fstream fid1("p.txt", std::ios::out);
    std::fstream fid2("p_l.txt", std::ios::out);
    std::fstream fid3("p_q.txt", std::ios::out);
    std::fstream fid4("p_lq.txt", std::ios::out);
    for (int i=0; i<M; ++i) {
        if (i > 0) fid1 << "\n";
        if (i > 0) fid2 << "\n";
        if (i > 0) fid3 << "\n";
        if (i > 0) fid4 << "\n";
        for (int j=0; j<N; ++j) {
            if (j > 0) fid1 << " ";
            if (j > 0) fid2 << " ";
            if (j > 0) fid3 << " ";
            if (j > 0) fid4 << " ";
            fid1 << phi[j + i*N];
            fid2 << phi_l[j + i*N];
            fid3 << phi_q[j + i*N];
            fid4 << phi_lq[j + i*N];
        }
    }
    delete [] phi;
    delete [] phi_l;
    delete [] phi_q;
    delete [] phi_lq;
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

