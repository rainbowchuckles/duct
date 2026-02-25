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
    if (argc < 3) {
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
read_config_file(argv + 1, gam, t, dt, m, u1, p1, Tv1, conv);

cout << "Inlet conditions: " << endl;
cout << endl;
cout << "p: " << " " << p1 << endl;
cout << "u: " << " " << u1 << endl;
cout << "Tv:" << " " << Tv1 << endl;
cout << endl;
cout << "---------------------" << endl;
cout << endl;

// setup the grid and area function
vector<double> x(m);
double dx = 0.001;
vector<double> l(N,0.0);
vector<double> Ai(N,1.0f);
vector<double> A(m,0.0f);
double dA;
double d;
int s;
bool flag;

// the area function
read_nozzle_file(argv + 3, l, Ai, s);

// create the grid in non-dimensional coordinates
for (int i = 1; i < m; i++){
        x[i] = x[i-1] + dx;
        dx = 1.00*dx;
}
// convert to dimensional
for (int i = 1; i < m; i++){
        x[i] *= l[s]/x[m-1];
	A[i] = 1.0 + 1.0*4.64*static_cast<float>(i)/static_cast<float>(m);
}
A[0] = 1.0;
// interpolate to find area for given x grid 
//linterp( l.data(), Ai.data(), s, x.data(), A.data(), m);

// read inlet file
vector<vector<double>> inlet;
vector<vector<double>> data;

read_inlet_file("inlet.dat",inlet);

// interpolate inlet file to match time grid
double Si[T][NSPMAX] = {0.0};
double y2[T];

std::vector<double> col3(T);
col3[0] = 0.0;
for (int k = 1; k < t; k++){col3[k] = col3[k-1] + dt;}

for (int o = 0; o < Ps+1; o++){
	vector<double> col0 = getColumn(inlet, 1  );
	vector<double> col1 = getColumn(inlet, o+2);
	linterp(col0.data(),col1.data(),inlet.size(),col3.data(),y2,t);
	for (int i = 0; i < t; ++i){
		Si[i][o]  = y2[i];
	if (o == Us){
		Si[i][o]  = -y2[i];
	}
	}
}

// e- N O N2 NO O2 
double rho = 0.0; 
vector<float> rho1(nsp,0.0);
rho1[3] = 0.00610;
rho1[5] = 0.00185;


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
double  F[NSPMAX]  = {0.0};
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
for (int i=0; i<m; i++){
for (int j=0; j<Ps+1; j++){S[k][i][j] = qp[j];}
}}

// the main loop
for (int k=0; k<t+1; k++){

//// apply the inlet boundary condition
//for (int j=0; j<Ps+1; j++){S[k][0][j] = Si[0][j];}
//
//if (k > 10000){
//for (int j=0; j<Ps+1; j++){S[k][0][j] = Si[k-10000][j];}
//}

for (int i=1; i<m; i++){

// 3. Evaluate G(S) = [0, ..., 0, 0, 0, p.dA/dx] -> geometric source term
flux1(  S[k][i],F,&nsp,wsp,asp,hfsp,rsp);

dA = A[i] - A[i-1];
dx = x[i] - x[i-1];

for (int j=0; j<X+1; j++){
	G[j] = 0.0;
	G[j] = -F[j]*dA/(dx*A[i]);
}
G[X] += S[k][i][Ps]*dA/(A[i]*dx);

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

for (int n = 0; n < X+1; n++){
	R[n] = F2[n] - F1[n];
	R[n] /= -dx;
	//R[n] += A[i]*Q[n];
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
if (k % 250 == 0){
	cout << k << " " << static_cast<float>(k)*dt << endl;
	write_cl(S, m, x, k, dt, nsp, wsp, "../out/cl.dat");
	write_outlet(S, m, x, t, dt, nsp, wsp, "../out/outlet.dat");
}
}



// output the centerline profile to text

write_cl(S, m, x, t-2, dt, nsp, wsp, "../out/cl.dat");
write_outlet(S, m, x, t, dt, nsp, wsp, "../out/outlet.dat");

cout << "Ran to completion. " << endl;


return 0;
}


