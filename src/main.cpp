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
#include "src.h"

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
for (int i = 0; i < m; i++){
	x[i] = i/static_cast<float>(m);
	x[i] = x[i]*l;
}

// the area function
vector<float> A(m,1.0f);

//for (int i =0; i<m; i++){
//	A[i] += 1*4.64*static_cast<float>(i)/m;
//}


// e- N O N2 NO O2  
vector<float> rho1(nsp,0e-3);
rho1[5] = 0.00687;
rho1[7] = 0.002086;

// the flow vector
// q(N,M,T)
// N: the primitive variables (rho, u, p)
// M: the streamwise points 
// T: the time steps
// this guy is defined in class.h

FlowField q(nv,m,t);

// initialise the flow field
for (int i =0; i<m; i++){
	for (int k= 0; k<t; k++){
		for (int j=0; j<nsp; j++){
		q(j,i,k) = rho1[j];
		}
		q(V,i,k) =  Tv1;
		q(U,i,k) =   u1;
		q(P,i,k) =   p1;
	}
}


// Euler updates
float dx;
float dA;
float du;
float dp;
float dTv;
float ye;
float rho = 0;
vector<float> drho(nsp);
int w;
float scale;

float pe = 0.0;  // the electron pressure
float cvv = 1.0; // the vibrational specific heat at constant volume

cout << "t, #  " << " " << "dt, s " << endl;
cout << endl;

// loop through time
for (int k = 0; k<t-1; k++){
	// print some indication of progress to the terminal
	if (k % 50 == 0 || k == 0 ){
		cout << setw(6) << k << " " 
		     << setw(6) << setprecision(3) << scientific << k*dt << endl;
	}
	// inlet condition
	for (int j=0; j<nsp; j++){
	q(j,0,k) = rho1[j];	
	}
	q(V,0,k) =   Tv1;	
	q(U,0,k) =   u1;	
	q(P,0,k) =   p1;                             
       	
 	// loop through space
	for (int i = 1; i<m; i++){
	dx   = x[i] - x[i-1];
	dA   = A[i] - A[i-1];
	rho = 0.0;
	for (int j=0; j<nsp; j++){
	 drho[j] = q(j,i,k) - q(j,i-1,k);
	 rho += q(j,i,k);
	}
	dTv  = q(V,i,k) - q(V,i-1,k);
	du   = q(U,i,k) - q(U,i-1,k);
	dp   = q(P,i,k) - q(P,i-1,k);
	
	// get the non-equilibrium source terms
	// only introduce them after the flow is steady
	//if (k*dt > l/u1){
	if (k > 1500){
	for (int o = 0; o < nsp; o++){qp[o]=q(o,i,k);}
	qp[V] = q(V,i,k);
	qp[U] = q(U,i,k);
	qp[P] = q(P,i,k);

	// continuity update
	for (int j=0; j<nsp; j++){
	q(j,i,k+1)  = dA/dx;         
	q(j,i,k+1) *= q(U,i,k)/A[i];
	q(j,i,k+1) += du/dx;        
	q(j,i,k+1) *= q(j,i,k);
	q(j,i,k+1) += q(U,i,k)*drho[j]/dx; 
	q(j,i,k+1) *= dt;

	q(j,i,k+1) = q(j,i,k) - q(j,i,k+1);
	}
	
	// vibrational energy update
	q(V,i,k+1)  = du/dx;
	q(V,i,k+1) += (q(U,i,k)/A[i])*(dA/dx);
	q(V,i,k+1) *= -pe;
	q(V,i,k+1) /= rho*cvv;
	q(V,i,k+1) += q(U,i,k)*dTv/dx;
	q(V,i,k+1) *= dt;

	q(V,i,k+1) = q(V,i,k) - q(V,i,k+1);

	// momentum update
	q(U,i,k+1)  = dp/dx;
	q(U,i,k+1) /= rho;      
	q(U,i,k+1) += q(U,i,k)*du/dx;
	q(U,i,k+1) *= dt;              
	
	q(U,i,k+1) = q(U,i,k) - q(U,i,k+1);              

	// energy update
	q(P,i,k+1) = dA/dx;         
	q(P,i,k+1) *= q(U,i,k)/A[i];
	q(P,i,k+1) += du/dx;        
	q(P,i,k+1) *= gam*q(P,i,k);   
	q(P,i,k+1) += q(U,i,k)*dp/dx;
	q(P,i,k+1) *= dt;              
	
	q(P,i,k+1) = q(P,i,k) - q(P,i,k+1);              
	
	src_c(&nsp,nelsp,ielsp,melsp,wsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr,qp,f);
	scale = static_cast<float>(k)/static_cast<float>(t);
	scale = min(static_cast<float>(1.0),scale);
	for (int o = 0; o < P+1; o++){q(o,i,k)= q(o,i,k)-f[o]*dt;}
	for (int p = 0; p < nsp; p++){ye += f[p];}
	
	}
	// check if steady state has been reached, exit if so
	// decide based on the exit Tv
	//if (k % 250 == 0 && k*dt > 1.1*l/u1 && k < t-1){
	if (k % 250 == 0 && k > 1600 && k < t-1){
		w = k;
		write_cl( q, m, x, w, nsp, wsp, "../out/cl.dat");
	//	if (abs(q(P,m-1,k) - q(P,m-1,k-100)) < conv){	
	//		w = k;
	//		goto post;} 
	}

	w = k;
	}

}

// output the centerline profile to text
post:
write_cl( q, m, x, w, nsp, wsp, "../out/cl.dat");
cout << endl;
cout << "Ran to completion. " << endl;


return 0;
}


