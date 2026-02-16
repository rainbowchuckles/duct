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
int T   = nsp+1;
int U   = nsp+2;

// read the user input file
read_input_file(argv + 1, gam, t, dt, m, l, u1, T1, Tv1, conv);

cout << "Inlet conditions: " << endl;
cout << endl;
cout << "T: " << " " << T1 << endl;
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
vector<float> rho1(nsp,0e-3);
rho1[3] = 0.00687;
rho1[5] = 0.002086;


// assemble the initial vector of primitive variables

for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
	for (int j=0; j<nsp; j++){qp[j] = rho1[j];}
	qp[V] = Tv1; 
	qp[U] = u1; 
	qp[T] = T1;
}}

// assemble the initial vector of conserved variables
//s2c_c(qp,qc,&nsp,wsp,asp,hfsp,rsp);

// the flow vector
// Q: N x M x T
// N: the primitive variables (rho_s,rhoev,rhou,rhoE)
// M: the streamwise points 
// T: the time steps

double Q[t+1][m+1][nv] = {0.0};

// the flux vectors

double F[nv]  = {0.0};

// initialise the solver
for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
for (int j=0; j<X+1; j++){Q[k][i][j] = qp[j]*A[i];}
}}

// the main loop
for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
// Recipe
// Given knowledge of the state vector S = [rhoi, ..., rhos, Tve, Ttr, u]
// 0. Calculate the static pressure from Boyle's law

boyle(p, S, &nsp, wsp, asp, hfsp, rsp);

// 1. Evaluate U(S) = [rhoi, ..., rhos, rhoev,rhoE,rhou] -> vector of conserved variables

s2u_c(S, U, Aus, &nsp, wsp, asp, hfsp, rsp);

// 2. Evaluate F(S) = [rhoiu, ..., rhosu, rhoevu,rhoHu,rhou^2 + p] -> flux of conserved variables

s2f_c(S, F, Afs, &nsp, wsp, asp, hfsp, rsp);

// 3. Evaluate G(S) = [0, ..., 0, 0, 0, p.dA/dx] -> geometric source term

G[X] = p*dA/dx;

// 4. Evaluate Q - > vector of thermochemical source terms
// 5. Aus = dU/dS, Afs = dF/dS - > Jacobians

src_c(Q);

// 6. For each streamwise cell i evaluate the numerical fluxes:
// F_{i+1/2} = 0.5 * ( F(S_L) + F(S_R) )
//             - 0.5 * alpha_{i+1/2} * ( U(S_R) - U(S_L) )
//
// where:
//   S_L  = state to the left of the interface
//   S_R  = state to the right of the interface
//   F(.) = physical flux function
//   U(.) = conserved variable vector
//   alpha_{i+1/2} = max wave speed at the interface

flux(F1, F2);

// 7. Form the residual in conserved variables
// R_{i,k} = -(1/dx)*(F_{i+1/2} - F_{i-1/2}) + Q

for (int n = 0; n < X+1; n++){
	R[n] = F2[n] - F1[n];
	R[n] /= -dx;
	R[n] += Q[n];
}

// 8. Solve the linear system for dS/dt
// dS/dt = A^-1 .R

// 9. Explicitly advance the state vector
// S_{i,k+1} = S_{i,k} + dt*dS/dt 
}}




// output the centerline profile to text
//post:
//write_cl( q, m, x, w, nsp, wsp, "../out/cl.dat");
//write_cl_mole( q, m, x, w, nsp, wsp, "../out/cl_mole.dat");
//cout << endl;
//cout << "Ran to completion. " << endl;


return 0;
}


