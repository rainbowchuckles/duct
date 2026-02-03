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


#include <fstream>
#include <string>
#include <stdexcept>

void read_input_file(char** filename,
                     float& gam,
                     int&   t,
                     float& dt,
                     int&   m,
                     float& l,
                     float& u1,
                     float& p1,
		     float& Tv1)
{
    // filename points to a C-string (e.g. argv + 1), so dereference it:
    const char* fname = *filename;      // same as filename[0]
    std::ifstream file(fname);
    if (!file) {
        throw std::runtime_error(std::string("Could not open input file: ") + fname);
    }

    float value;
    std::string name;

    while (file >> value >> name) {
        if (name == "gam") {
            gam = value;
        } else if (name == "t") {
            t = static_cast<int>(value);
        } else if (name == "dt") {
            dt = value;
        } else if (name == "m") {
            m = static_cast<int>(value);
        } else if (name == "l") {
            l = value;
        } else if (name == "u1") {
            u1 = value;
        } else if (name == "p1") {
            p1 = value;
        } else if (name == "Tv1") {
            Tv1 = value;
        }
    }
}

