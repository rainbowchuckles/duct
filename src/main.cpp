#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

// Park vars

#define NELMAX 8
#define NSPMAX 64
#define NTRMAX 4
#define NACOEF 12
#define NRNMAX 128

#include "class.h"
#include "gas.h"
#include "post.h"

extern "C" void tcw_c(
    int *nsp,

    double nelsp[NSPMAX],

    // Fortran: ielsp(NELMAX, NSPMAX)
    double ielsp[NSPMAX][NELMAX],

    // Fortran: melsp(NELMAX, NSPMAX)
    double melsp[NSPMAX][NELMAX],

    double wsp[NSPMAX],

    // Fortran: rsp(4, NSPMAX)
    double rsp[NSPMAX][4],

    // Fortran: asp(NACOEF, NTRMAX, NSPMAX)
    double asp[NSPMAX][NTRMAX][NACOEF],

    double hfsp[NSPMAX],

    // Fortran: mw(3, NSPMAX, NSPMAX)
    double mw[NSPMAX][NSPMAX][3],

    // Fortran: cs(2, NSPMAX, NSPMAX)
    double cs[NSPMAX][NSPMAX][2],

    double *diss,
    double *inz,

    // Fortran: apb(3, NSPMAX)
    double apb[NSPMAX][3],

    int nrn[2],

    // Fortran: nsprn(4, NRNMAX)
    int nsprn[NRNMAX][4],

    // Fortran: isprn(NSPMAX, NRNMAX)
    int isprn[NRNMAX][NSPMAX],

    // Fortran: msprn(NSPMAX, NRNMAX)
    int msprn[NRNMAX][NSPMAX],

    // Fortran: ktbrn(NRNMAX)
    int ktbrn[NRNMAX],

    // Fortran: xtbrn(NSPMAX, NRNMAX)
    int xtbrn[NRNMAX][NSPMAX],

    // Fortran: arr(64, NSPMAX)
    double arr[NSPMAX][64]
);

extern "C" void src_c(
    int *nsp,

    double nelsp[NSPMAX],

    double ielsp[NSPMAX][NELMAX],
    double melsp[NSPMAX][NELMAX],

    double wsp[NSPMAX],

    double rsp[NSPMAX][4],
    double asp[NSPMAX][NTRMAX][NACOEF],

    double hfsp[NSPMAX],

    double mw[NSPMAX][NSPMAX][3],
    double cs[NSPMAX][NSPMAX][2],

    double *diss,
    double *inz,

    double apb[NSPMAX][3],

    int nrn[2],
    int nsprn[NRNMAX][4],
    int isprn[NRNMAX][NSPMAX],
    int msprn[NRNMAX][NSPMAX],

    int ktbrn[NRNMAX],
    int xtbrn[NRNMAX][NSPMAX],

    double arr[NSPMAX][64],

    double *q,
    double *f
);


using namespace std;

// the main program

int main() {

double nelsp[NSPMAX];
double ielsp[NSPMAX][NELMAX];
double melsp[NSPMAX][NELMAX];
double wsp[NSPMAX]; // Fix
double rsp[NSPMAX][4];
double asp[NSPMAX][NTRMAX][NACOEF];
double hfsp[NSPMAX];
double mw[NSPMAX][NSPMAX][3];
double cs[NSPMAX][NSPMAX][2];
double diss[NSPMAX];
double inz[NSPMAX];
double apb[NSPMAX][3];
int    nrn[2];
int    nsprn[NRNMAX][4];
int    isprn[NRNMAX][NSPMAX];
int    msprn[NRNMAX][NSPMAX];
int    ktbrn[NRNMAX];
int    xtbrn[NRNMAX][NSPMAX];
double arr[NSPMAX][64];

// species
int nsp = 1;

// bs for src
double qp[NSPMAX];
double f[NSPMAX];

// Retrieve thermo stuff
tcw_c(&nsp,nelsp,ielsp,melsp,wsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr);

// species
int nv  = nsp + 3;
int V = nsp;
int U = nsp+1;
int P = nsp+2;

// gas
// the ratio of specific heats
float gam = 1.4;

// time       
// t is the number of time steps
int t = 3500;

// dt is the size of the time step
float dt = 1.0e-6;

// grid
// m is the number of grid points
int m = 400;

// l is the length of the domain
float l = 2.2;

// x is the vector of physical grid locations
vector<float> x(m);

// populate the grid vector
	for (int i = 0; i < m; i++){
		x[i] = i/static_cast<float>(m);
		x[i] = x[i]*l;
	}

// the area function
vector<float> A(m,1.0f);

//	for (int i =0; i<m; i++){
//		A[i] = A[i] + 4.64*static_cast<float>(i)/m;
//			}


// the inlet boundary condition
vector<float> rho1(nsp,0e-3);
float u1   = 3500;
float p1   = 7000;
float Tv1  = 2500;

// e- N O N2 NO O2  
rho1[3] = 1e-3;
rho1[5] = 1e-3;

// the flow vector
// q(N,M,T)
// N: the primitive variables (rho, u, p)
// M: the streamwise points 
// T: the time steps

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
float rho = 0;
vector<float> drho(nsp);

float pe = 0.0;  // the electron pressure
float cvv = 200.0; // the vibrational specific heat at constant volume

	// loop through time
	for (int k = 0; k<t-1; k++){
		// print some indication of progress to the terminal
		if (k % 50 == 0){cout << k << endl;}
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
		
		// get the non-equilibrium source terms
		// only introduce them after the flow is steady
		if (k > 1500){
		for (int o = 0; o < nsp; o++){qp[o]=q(o,i,k);}
		qp[V] = q(V,i,k);
		qp[U] = q(U,i,k);
		qp[P] = q(P,i,k);
		
		src_c(&nsp,nelsp,ielsp,melsp,wsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr,qp,f);

		for (int o = 0; o < P+1; o++){q(o,i,k)= q(o,i,k)-f[o]*dt;}
		}
		}
	}
	write_cl( q, m, x, t, nsp, wsp, "../out/cl.dat");


  	return 0;
}


