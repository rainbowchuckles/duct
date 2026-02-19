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

#include "macro.h"

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
double dx = 0.001;

// populate the grid vector
for (int i = 1; i < m; i++){
	x[i] = x[i-1] + dx;
	dx = 1.02*dx;
}
for (int i = 0; i < m; i++){x[i] *= l/x[m-1];}

// the area function
vector<float> A(m,1.0f);

for (int i =0; i<m; i++){
	A[i] += 1*4.64*static_cast<float>(i)/m;
}


// e- N O N2 NO O2 
double rho = 0.0; 
vector<float> rho1(nsp,0.0);
rho1[5] = 0.00687;
rho1[7] = 0.002086;


// assemble the initial vector of primitive variables

for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
	for (int j=0; j<nsp; j++){qp[j] = rho1[j];}
	qp[Vs] = Tv1; 
	qp[Us] = u1; 
	qp[Ps] = p1;
}}

// the flow vector
double S[T][M][NSPMAX] = {0.0};
double U[T][M][NSPMAX] = {0.0};

// the flux vectors
double F1[NSPMAX]  = {0.0};
double F2[NSPMAX] = {0.0};

// source vectors
double G[NSPMAX]  = {0.0};
double Q[NSPMAX]  = {0.0};

// residual vector
double R[NSPMAX]  = {0.0};
double C[NSPMAX]  = {0.0};

// initialise the solver
for (int k=0; k<t+1; k++){
for (int i=0; i<m+1; i++){
for (int j=0; j<nsp; j++){S[k][i][j] = qp[j];}
S[k][i][Vs] = qp[Vs];
S[k][i][Us] = qp[Us];
S[k][i][Ps] = qp[Ps];
}}
for (int k=0; k<t+1; k++){
for (int i=1; i<m+1; i++){
S[k][i][Ps] = 24627;
}}

// the main loop
for (int k=0; k<t+1; k++){
for (int i=1; i<m; i++){

// 3. Evaluate G(S) = [0, ..., 0, 0, 0, p.dA/dx] -> geometric source term
double dA = A[i+1] - A[i];
G[X] = S[k][i][Ps]*dA/dx;

// 4. Evaluate Q - > vector of thermochemical source terms
for (int j=0; j<Ps+1; j++){qp[j] = S[k][i][j];}

src_c(&nsp,nelsp,ielsp,melsp,wsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr,qp,Q,J);

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

flux(S[k][i-1],  S[k][i],F1,&nsp,wsp,asp,hfsp,rsp);
flux(  S[k][i],S[k][i+1],F2,&nsp,wsp,asp,hfsp,rsp);

// supersonic outlet
if (i == m-1){
for (int n = 0; n<X+1; n++){
		F2[n] = F1[n];
		}
}

// 7. Form the residual in conserved variables
// R_{i} = -(1/dx)*(F_{i+1/2} - F_{i-1/2}) + Q

dx = x[i+1] - x[i];
for (int n = 0; n < X+1; n++){
	R[n] =  F2[n] - F1[n];
	R[n] /= -dx;
	R[n] += Q[n];
	R[n] += G[n];
}

// 8. Solve the linear system for dS/dt
// dS/dt = J^-1 .R
// J obtained from source is apc (i.e. exactly the matrix we want) 

// Matrix multiplication: C = J * R

for (int j = 0; j < Ps+1; j++) {
	C[j] = 0.0;
 for (int n=0; n<X+1; n++){
            C[j] += J[n][j] * R[n];
	}
}

// 9. Explicitly advance the state vector
// S_{i,k+1} = S_{i,k} + dt*dS/dt 
for (int n=0; n<Ps+1; n++){
	S[k+1][i][n] = S[k][i][n] + dt*C[n];
}


// supersonic outlet
if (i == m-1){
for (int n = 0; n<X+1; n++){
		S[k+1][i][n] = S[k+1][i-1][n];
		}
}

}
if (k % 100 == 0){
	cout << k << " " << static_cast<float>(k)*dt << endl;
	write_cl(S, m, x, k, nsp, wsp, "../out/cl.dat");
}
}



// output the centerline profile to text
//post:
write_cl(S, m, x, t-2, nsp, wsp, "../out/cl.dat");
//write_cl_mole( q, m, x, w, nsp, wsp, "../out/cl_mole.dat");
//cout << endl;
//cout << "Ran to completion. " << endl;


return 0;
}


