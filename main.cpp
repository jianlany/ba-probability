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
        "\t\t--q_width_factor: same but on angle direction. Default: 1.5\n"
        "\t\t--d_range: the upper and lower limit of bond. Default: 2.0  2.9\n"
        "\t\t--theta_range: the upper and lower limit of angle. Default: 0  180\n"
        "\t\t--output-file: the path of the output badf file. Default: p.txt\n"
        "\t\t--help: print this message.\n";
struct BA {
    double d1, d2;
    double q;
};
using BA_data = std::vector<BA>;

BA_data read_BA_file(std::string path, double wq);


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Usage:\n  ./ba-probability [path] (options)\n.";
        return 1;
    }
    std::string ba_path = argv[1];
    std::string output_path = "p.txt";

    int nl=100, nq=100;
    auto xlo = 2.0,  xhi = 2.9;
    auto ylo = 0.0, yhi = 180.0;
    // How many grid spacings the distribution is spread over.
    double d_width_factor = 1.5;
    double q_width_factor = 1.5;
    // Loop over remaining arguments and set optional parameters.
    bool entropy_factor = false;
    for (int i=2; i<argc; ++i) {
        if (!strcmp(argv[i], "--help")) {
            std::cout << help_msg;
            exit(0);
        }
        else if (!strcmp(argv[i], "--num_theta") && ++i < argc) {
            nq = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--num_length") && ++i < argc) {
            nl = atoi(argv[i]);
        }
        else if (!strcmp(argv[i], "--d_width_factor") && ++i < argc) {
            d_width_factor = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--q_width_factor") && ++i < argc) {
            q_width_factor = atof(argv[i]);
        }
        else if (!strcmp(argv[i], "--d_range") && (i+2) < argc) {
            xlo = atof(argv[++i]);
            xhi = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--theta_range") && (i+2) < argc) {
            ylo = atof(argv[++i]);
            yhi = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--output-file") && ++i < argc) {
            output_path = argv[i];
        }
        else if (!strcmp(argv[i], "--entropy")) {
            entropy_factor = true;
        }
        else{
            std::cerr << "Invalid input arguments: " << argv[i] << " \n";
            std::cerr << help_msg;
            exit(1);
        }
    }
    // Compute grid spacing.
    const auto dx = (xhi-xlo) / (nl-1);
    const auto dy = (yhi-ylo) / (nq-1);
    const auto N = nq * nl;
    // kernel width is 1.5 grid spacings.
    const auto wx = d_width_factor*dx;
    const auto wy = q_width_factor*dy;
    auto ba_data = read_BA_file(ba_path, wy);
    const double *data = &ba_data[0].d1;
    int n_ba = ba_data.size();
    const auto Z = 1.0 / (2.0 * pi * wx*wy) / (2.0 * ba_data.size());
    double *p = new double[N];
    double *p_l = new double[N];
    double *p_q = new double[N];
    double *p_lq = new double[N];
    for (int i=0; i<N; ++i) {
        p[i] = p_l[i] = p_q[i] = p_lq[i] = 0.0;
    }
    #pragma acc data copy(data[:3*n_ba], p[:N], p_l[:N], p_q[:N], p_lq[:N], entropy_factor)
    {
    double entropy_coeff = 1.0;
    #pragma acc parallel loop collapse(2)
    for (int j=0; j<nl; ++j) {
        for (int k=0; k<nq; ++k) {
            for (int i=0; i<n_ba; ++i) {
                auto d1 = data[3*i];
                auto d2 = data[3*i+1];
                auto q = data[3*i+2];
                auto xj = xlo + dx*j;
                auto yk = ylo + dy*k;
                auto x1 = (xj-d1) / wx;
                auto x2 = (xj-d2) / wx;
                auto y = (yk-q) / wy;
                if (entropy_factor) entropy_coeff = 1.0/fabs(sin(q/180.0*pi));
                auto exp1 = exp(-0.5*(x1*x1 + y*y))*entropy_coeff ;
                auto exp2 = exp(-0.5*(x2*x2 + y*y))*entropy_coeff ;

                p[k + j*nq] += (Z*exp1 + Z*exp2);
                p_l[k + j*nq] += Z*exp1*(-x1/wx) + Z*exp2*(-x2/wx);
                p_q[k + j*nq] += Z*exp1*(-y/wy) + Z*exp2*(-y/wy);
                p_lq[k + j*nq] += Z*exp1*(y/wy)*(x1/wx) + Z*exp2*(y/wy)*(x2/wx);
            }
        }
    }
    }

    auto sum = 0.0;
    for (int i=0; i<nl-1; ++i) {
        for (int j=0; j<nq-1; ++j) {
            auto p1 = p[j   + i*nq];
            auto p2 = p[j   + (i+1)*nq];
            auto p3 = p[j+1 + (i+1)*nq];
            auto p4 = p[j+1 + i*nq];
            sum += 0.25*(p1+p2+p3+p4) * dx * dy;
        }
    }
    std::cout << "Sum = " << sum << "\n";

    int j = nq-1;
    // Set the value at 180 same as the value at 180 - dq
    for (int i=0; i<nl; ++i){
         p[j + i*nq]    = p[j + i*nq - 1];
         p_l[j + i*nq]  = p_l[j + i*nq - 1];
         p_q[j + i*nq]  = 0;
         p_lq[j + i*nq] = 0;
    }

    // Renormalized the distribution
    for (int i=0; i<nl; ++i) {
        for (int j=0; j<nq; ++j) {
             p[j   + i*nq]    /= sum;
             p_l[j   + i*nq]  /= sum;
             p_q[j   + i*nq]  /= sum;
             p_lq[j   + i*nq] /= sum;
        }
    }

    // Sum should be close to 1.0.
    std::fstream fid(output_path, std::ios::out);

    fid << "# Probability: " << ba_path << "\n";
    fid << "# Grid size: " << nl << " " << nq << "\n";
    fid << "# (bond length) xlo xhi: " << xlo << " " << xhi << "\n";
    fid << "# (bond angle)  ylo yhi: " << ylo << " " << yhi << "\n";
    auto write_table = [&](const double *x, const char *label) {
        fid << "# " << label << "\n";
        fid.precision(16);
        for (int j=0; j<nq; ++j) {
            for (int i=0; i<nl; ++i) {
                if (i > 0) fid << " ";
                fid << x[j + i*nq];
            }
            fid << "\n";
        }
    };
    write_table(p, "p");
    write_table(p_l, "dp/dl");
    write_table(p_q, "dp/dq");
    write_table(p_lq, "dp^2/dldq");
    delete [] p;
    delete [] p_l;
    delete [] p_q;
    delete [] p_lq;
}


// Reads the BA datafile output from CG-distributions.  Each row of the file 
// should contain: [bond_len1, bond_len2, bond_angle].
BA_data read_BA_file(std::string path, double wq) {
    BA_data ba_data;
    std::ifstream fid(path, std::ios::in);
    if (!fid) {
        std::cerr << "Could not open " << path << ".\n";
        exit(1);
    }
    while (fid) {
        BA ba;
        fid >> ba.d1 >> ba.d2 >> ba.q;
        if (ba.q > 179.99) continue;
        if (fid) {
            ba_data.push_back(ba);
        }
        // If there the data point is close to 180 degree
        // Add an extra one on the other side to avoid fall-off at 180 degree
        // Not necessary at 0 degree since it's not physically possible.
        if ((180.0 - ba.q) < 5*wq){
            ba.q = 360 - ba.q;
            ba_data.push_back(ba);
        }
    }
    return ba_data;
}

