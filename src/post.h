#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;


void write_flowfield(
    const FlowField& q,
    int m,
    int t,
    const std::string& filename
) {
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << '\n';
        return;
    }

    // header (commented for plotting tools)
    file << "# i  k  rho  u  p\n";

    file << std::scientific << std::setprecision(10);

    for (int k = 0; k < t; ++k) {
        for (int i = 0; i < m; ++i) {
            file
                << i << " "
                << k << " "
                << q(0,i,k) << " "   // rho
                << q(1,i,k) << " "   // u
                << q(2,i,k) << "\n"; // p
        }
        file << "\n";  // blank line between time slices (optional)
    }

    file.close();
}

void write_outlet(
    const FlowField& q,
    int m,
    int t,
    const std::string& filename
) {
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << '\n';
        return;
    }

    // header (commented for plotting tools)
    file << "# i  k  rho  u  p\n";

    file << std::scientific << std::setprecision(10);

    for (int k = 0; k < t; ++k) {
            file
                << m-1 << " "
                << k << " "
                << q(0,m-1,k) << " "   // rho
                << q(1,m-1,k) << " "   // u
                << q(2,m-1,k) << "\n"; // p
    }
    file.close();
}

void write_cl(
    const FlowField& q,
    int m,
    int t,
    const std::string& filename
) {
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << '\n';
        return;
    }

    // header (commented for plotting tools)
    file << "# i  k  rho  u  p\n";

    file << std::scientific << std::setprecision(10);

    for (int i = 0; i < m; ++i) {
        file
            << i << " "
            << t-1 << " "
            << q(0,i,t-1) << " "   // rho
            << q(1,i,t-1) << " "   // u
            << q(2,i,t-1) << "\n"; // p
    }

    file.close();
}
