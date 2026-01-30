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
    int nsp,
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
                << k << " "; 
		for (int j =0; j<nsp; j++){
                file << q(j,i,k) << " ";  // rho
		}
            file
		<< q(nsp,i,k) << " "   // u
		<< q(nsp+1,i,k) << " "        
                << q(nsp+2,i,k) << "\n"; // p
        }
        file << "\n";  // blank line between time slices (optional)
    }

    file.close();
}

void write_outlet(
    const FlowField& q,
    int m,
    int t,
    int nsp,
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
                << k << " "; 
		for (int j =0; j<nsp; j++){
                file << q(j,m-1,k) << " ";  // rho
		}
            file
		<< q(nsp,m-1,k) << " "   // u
		<< q(nsp+1,m-1,k) << " "   // u
                << q(nsp+2,m-1,k) << "\n"; // p
    }
    file.close();
}

void write_cl(
    const FlowField& q,
    int m,
    vector<float> x,
    int t,
    int nsp,
    double wsp[NSPMAX],
    const std::string& filename
) {
    float tt;
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error opening file: " << filename << '\n';
        return;
    }

    // header (commented for plotting tools)
    file << "# i  k  rho  u  p t\n";

    file << std::scientific << std::setprecision(10);

    for (int i = 0; i < m; ++i) {
            file
                << x[i] << " "
                << t-1 << " "; 
		for (int j =0; j<nsp; j++){
                file << q(j,i,t-1) << " ";  // rho
		}
            file
		<< q(nsp,i,t-1) << " "   // u
		<< q(nsp+1,i,t-1) << " "   // u
                << q(nsp+2,i,t-1) << " "; // p
            tt = calt(q,i,t-1,nsp,wsp);
            file
		<< tt << "\n";   // u
    }

    file.close();
}
