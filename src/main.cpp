// external libraries 
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include <cstring>

// variable definitions
#include "size.h"
#include "class.h"
#include "gas.h"
#include "io.h"
#include "func.h"

using namespace std;

// the main program

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input_filename\n";
        return 1;
    }

#include "vars.h"

// initialise the non-equilibrium model
// tcw_c is a wrapper for the OCEAN Fortran functions in tcw.F
const char* inp = argv[2];
int length = strlen(inp);

tcw_c(inp,&nsp,nelsp,ielsp,melsp,wsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr);

// nv is the number of variables
// nsp species + Tv + u + p 
int nv  = nsp + 4;
int V   = nsp;
int U   = nsp+1;
int P   = nsp+2;

int I   = nsp;
int E   = nsp+1;
int X   = nsp+3;

// read the user input file
read_input_file(argv + 1, gam, t, dt, m, l, u1, p1, Tv1, conv);

cout << "Inlet conditions: " << endl;
cout << endl;
cout << "p: " << " " << p1 << endl;
cout << "u: " << " " << u1 << endl;
cout << "Tv:" << " " << Tv1 << endl;
cout << endl;
cout << "---------------------" << endl;
cout << endl;

// x is the vector of physical grid locations
vector<float> x(m);

// populate the grid vector
for (int i = 0; i < m; i++){x[i] = static_cast<float>(i)/static_cast<float>(m);}
for (int i = 0; i < m; i++){x[i] *= l;}

// the area function
vector<float> A(m,1.0f);

//for (int i =0; i<m; i++){
//	A[i] += 1*4.64*static_cast<float>(i)/m;
//}


// e- N O N2 NO O2 
double rho = 0.0; 
double p   = 0.0; 
vector<float> rho1(nsp,0e-3);
rho1[3] = 0.00687;
rho1[5] = 0.002086;


// assemble the initial vector of primitive variables

for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
	for (int j=0; j<nsp; j++){qp[j] = rho1[j];}
	qp[V] = Tv1; 
	qp[U] = u1; 
	qp[P] = p1; 

}}

// assemble the initial vector of conserved variables
s2c_c(qp,qc,&nsp,wsp,asp,hfsp,rsp);

// the flow vector
// Q: N x M x T
// N: the primitive variables (rho_s,rhoev,rhou,rhoE)
// M: the streamwise points 
// T: the time steps

double Q[t][m][NSPMAX] = {0.0};

// the left and right flux vectors

double Fl[NSPMAX] = {0.0};
double Fr[NSPMAX] = {0.0};


// initialise the solver
for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
for (int j=0; j<X+1; j++){Q[k][i][j] = qc[j]*A[i];}
}}


// the main loop
for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
	// convert to primitive for flux calculation	
	// densities
	rho = 0.0;
	for (int j=0; j<nsp; j++){
		qp[j] = Q[k][i][j]/A[i]; 
		rho += qp[j];}
	// vibrational energy
	qp[V] = Q[k][i][E]/rho;
	// velocity
	qp[U] = Q[k][i][X]/rho;
        // total enthalpy (H) from internal energy (E)
	qp[I] = Q[k][i][I]/rho;
	cout << qp[I] << endl;
	p = 0.4*rho*(qp[I] - 0.5*qp[U]*qp[U]);
	qp[I] += p/rho;
	cout << p   << endl;
	cout << qp[I] << endl;
	exit(0);
	// form the left flux vector F_L
	
}}




// output the centerline profile to text
//post:
//write_cl( q, m, x, w, nsp, wsp, "../out/cl.dat");
//write_cl_mole( q, m, x, w, nsp, wsp, "../out/cl_mole.dat");
//cout << endl;
//cout << "Ran to completion. " << endl;


return 0;
}


