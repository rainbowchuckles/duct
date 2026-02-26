#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;

void write_cl(
    double S[T][M][NSPMAX], 
    int m,
    vector<double> x,
    int t,
    float dt,
    int nsp,
    double wsp[NSPMAX],
    const std::string& filename 
) {
    float tt;
    float rho;
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
                << static_cast<float>(t-1)*dt << " "; 
		for (int j =0; j<nsp; j++){
                file << S[0][i][j] << " ";  // rho
		}
            file
		<< S[0][i][nsp] << " "   // u
		<< S[0][i][nsp+1] << " "   // u
                << S[0][i][nsp+2] << " "; // p
            tt = calt(S[0][i],nsp,wsp);
            file
		<< tt << " ";
	    rho = 0.0;
	    for (int o = 0; o <nsp; ++o){rho += S[0][i][o];}
	    file
		<< rho << "\n";
    }

    file.close();
}




// write temporal outlet profile

#include <filesystem>  // C++17

void write_outlet(
    double S[T][M][NSPMAX],
    int m,
    vector<double> x,
    int t,
    float dt,
    int nsp,
    double wsp[NSPMAX],
    const std::string& filename
) {
    float tt;
    float rho;

    bool file_exists = std::filesystem::exists(filename);

    std::ofstream file(filename, std::ios::app);

    if (!file) {
        std::cerr << "Error opening file: " << filename << '\n';
        return;
    }

    // Write header only if file didn't exist before
    if (!file_exists) {
        file << "# i  k  rho  u  p t\n";
    }

    file << std::scientific << std::setprecision(10);

    // write current step
    file
        << x[1] << " "
        << static_cast<float>(t) * dt << " ";

    for (int j = 0; j < nsp; j++) {
        file << S[0][m-1][j] << " ";
    }

    file
        << S[0][m-1][nsp]   << " "
        << S[0][m-1][nsp+1] << " "
        << S[0][m-1][nsp+2] << " ";

    tt = calt(S[0][m-1], nsp, wsp);
    file << tt << " ";

    rho = 0.0;
    for (int o = 0; o < nsp; ++o) {
        rho += S[0][m-1][o];
    }

    file << rho << "\n";

    file.close();
}



#include <fstream>
#include <string>
#include <stdexcept>

void read_config_file(char** filename,
                     float& gam,
                     int&   t,
                     float& dt,
                     int&   m,
                     float& u1,
                     float& p1,
		     float& Tv1,
		     float& conv)
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
        } else if (name == "u1") {
            u1 = value;
        } else if (name == "p1") {
            p1 = value;
        } else if (name == "Tv1") {
            Tv1 = value;
        } else if (name == "conv") {
            conv = value;
        }
    }
}

void read_thermo_files(char** filename,
		       string &chmf,
		       string &rcnf,
		       string &mw_file,
		       string &cs_file,
		       string &mst_file,
		       string &diss_file,
		       string &ion_file,
		       string &apb_file,
		       string &thrf,
		       string &colpth)
{

    // filename points to a C-string (e.g. argv + 1), so dereference it:
    const char* fname = *filename;      // same as filename[0]
    std::ifstream file(fname);
    if (!file) {
        throw std::runtime_error(std::string("Could not open input file: ") + fname);
    }

    std::string value;
    std::string name;

    while (file >> value >> name) {
        if (name == "chmf") {
            chmf = value.c_str();
        } else if (name == "rcnf") {
            rcnf = value;
        } else if (name == "mw_file") {
            mw_file = value;
        } else if (name == "cs_file") {
            cs_file = value;
        } else if (name == "mst_file") {
            mst_file = value;
        } else if (name == "diss_file") {
            diss_file = value;
        } else if (name == "ion_file") {
            ion_file = value;
        } else if (name == "apb_file") {
            apb_file = value;
        } else if (name == "thrf") {
            thrf = value;
        } else if (name == "colpth") {
            colpth = value;
        }
    }
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

bool read_inlet_file(const std::string& filename,
                   std::vector<std::vector<double>>& data)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    std::string line;

    // Skip header line
    std::getline(file, line);

    // Skip blank line (if present)
    std::getline(file, line);

    // Read numeric data
    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        std::stringstream ss(line);
        std::vector<double> row;
        double value;

        while (ss >> value)
        {
            row.push_back(value);
        }

        if (!row.empty())
            data.push_back(row);
    }

    file.close();
    return true;
}


// read nozzle geometry
void read_nozzle_file(char** filename,
                      vector<double>& x,
                      vector<double>& A,
                      int& s)
{
    // filename points to a C-string (e.g. argv + 1), so dereference it:
    const char* fname = *filename;      // same as filename[0]
    std::ifstream file(fname);
    if (!file) {
        throw std::runtime_error(std::string("Could not open geometry file: ") + fname);
    }

    std::vector<std::vector<double>> data(2);  // 2 rows
    double col1, col2;

    // read x and y points into data
    int i = 0;
    while (file >> col1 >> col2) {
            x[i] = col1;
            A[i] = 3.1415*col2*col2;
            i++;
    }
    s = i-1;
    file.close();

    return;

}
