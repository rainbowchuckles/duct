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
int nv  = nsp + 3;
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
double rhol = 0.0; 
double rhor = 0.0; 
double p   = 0.0; 
double pl   = 0.0; 
double pr   = 0.0; 
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

double Q[t+1][m+1][nv] = {0.0};

// the flux vectors

double F[nv]  = {0.0};
double Fl[nv] = {0.0};
double Fr[nv] = {0.0};

double S[nv] = {0.0};

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
	rhol = 0.0;
	for (int j=0; j<nsp; j++){
		qpl[j] = Q[k][i][j]/A[i]; 
		rhol += qpl[j];}
	for (i+1nt j=0; j<nsp; j++){
		qpl[j] = Q[k][i+1][j]/A[i+1]; 
		rhor += qpl[j];}
	// vibrational energy
	qpl[V] = Q[k][i][E]/rhol;
	qpr[V] = Q[k][i+1][E]/rhor;
	// velocity
	qpl[U] = Q[k][i][X]/rhol;
	qpr[U] = Q[k][i+1][X]/rhor;
        // total enthalpy (H) from internal energy (E)
	qpl[I] = Q[k][i][I]/rhol;
	pl = 0.4*rhol*(qpl[I] - 0.5*qpl[U]*qpl[U]);
	qpl[I] += pl/rhol;
	qpr[I] = Q[k][i+1][I]/rhor;
	pr = 0.4*rhor*(qpr[I] - 0.5*qpr[U]*qpr[U]);
	qpr[I] += pr/rhor;

	// form the left and right flux vectors
	// continuity	
	for (int j=0; j<nsp; j++){Fl[j] = qpl[j]*qpl[U]*(A[i+1]+A[i])/2.0:}
	for (int j=0; j<nsp; j++){Fr[j] = qpr[j]*qpr[U]*(A[i+1]+A[i])/2.0:}
	// vibrational energy
	Fl[E] = rhol*qpl[U]*qpl[V]*(A[i+1]+A[i])/2.0;
	Fr[E] = rhor*qpr[U]*qpr[V]*(A[i+1]+A[i])/2.0;
	// momentum
	Fl[X] = (rhol*qpl[U]*qpl[U] + pl)*(A[i+1]+A[i])/2.0;
	Fr[X] = (rhor*qpr[U]*qpr[U] + pr)*(A[i+1]+A[i])/2.0;
	// total energy
	Fl[I] = rhol*qpl[U]*qpl[I]*(A[i+1]+A[i])/2.0;
	Fr[I] = rhor*qpr[U]*qpr[I]*(A[i+1]+A[i])/2.0;

	// area variation source term
        S[X]  = (A[i+1]+A[i])/2.0; 
        S[X] -= (A[i]+A[i-1])/2.0; 
	S[X] *= pl;

	// discrete update
	
}}




// output the centerline profile to text
//post:
//write_cl( q, m, x, w, nsp, wsp, "../out/cl.dat");
//write_cl_mole( q, m, x, w, nsp, wsp, "../out/cl_mole.dat");
//cout << endl;
//cout << "Ran to completion. " << endl;


return 0;
}


