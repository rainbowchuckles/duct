#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

#include "class.h"
#include "gas.h"
#include "post.h"

// Park vars

#define NELMAX 8
#define NSPMAX 64
#define NTRMAX 4
#define NACOEF 12
#define NRNMAX 128

extern "C" void tcw_c(int *nsp, double *nelsp, double ielsp[][NSPMAX], double melsp[][NSPMAX], 
		      double *wsp,
		      double rsp[][NSPMAX], double asp[][NTRMAX][NSPMAX], double *hfsp, double mw[][NSPMAX][NSPMAX], 
		      double cs[][NSPMAX][NSPMAX],
		      double *diss, double *inz, double apb[][NSPMAX], int *nrn, int nsprn[][NRNMAX], 
		      int isprn[][NRNMAX], int msprn[][NRNMAX], int *ktbrn, int xtbrn[][NRNMAX], 
		      double arr[][NSPMAX]);

extern "C" void src_c(int *nsp, double *nelsp, double ielsp[][NSPMAX], double melsp[][NSPMAX], 
		      double *wsp,
		      double rsp[][NSPMAX], double asp[][NTRMAX][NSPMAX], double *hfsp, double mw[][NSPMAX][NSPMAX], 
		      double cs[][NSPMAX][NSPMAX],
		      double *diss, double *inz, double apb[][NSPMAX], int *nrn, int nsprn[][NRNMAX], 
		      int isprn[][NRNMAX], int msprn[][NRNMAX], int *ktbrn, int xtbrn[][NRNMAX], 
		      double arr[][NSPMAX], double *q, double *f);
using namespace std;

// the main program

int main() {


double nelsp[NSPMAX];
double ielsp[NELMAX][NSPMAX];
double melsp[NELMAX][NSPMAX];
double wwsp[NSPMAX]; // Fix
double rsp[4][NSPMAX];
double asp[NACOEF][NTRMAX][NSPMAX];
double hfsp[NSPMAX];
double mw[3][NSPMAX][NSPMAX];
double cs[2][NSPMAX][NSPMAX];
double diss[NSPMAX];
double inz[NSPMAX];
double apb[3][NSPMAX];
int    nrn[2];
int    nsprn[4][NRNMAX];
int    isprn[NSPMAX][NRNMAX];
int    msprn[NSPMAX][NRNMAX];
int    ktbrn[NSPMAX];
int    xtbrn[NSPMAX][NRNMAX];
double arr[64][NSPMAX];

// species
int nsp =12;
int nv  = nsp + 3;
int V = nsp;
int U = nsp+1;
int P = nsp+2;

// bs for src
double qp[NSPMAX];
double f[NSPMAX];

tcw_c(&nsp,nelsp,ielsp,melsp,wwsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr);
src_c(&nsp,nelsp,ielsp,melsp,wwsp,rsp,asp,hfsp,mw,cs,diss,inz,apb,nrn,nsprn,isprn,msprn,ktbrn,xtbrn,arr,qp,f);

for (int i=0; i<32; i++){cout << wwsp[i] << endl;}


exit(0);




// gas
// the ratio of specific heats
float gam = 1.4;

// species weights (to come from OCEAN later)
vector<float> wsp(nsp,1.e-3);

// time       
// t is the number of time steps
int t = 15000;

// dt is the size of the time step
float dt = 1.e-7;

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

	for (int i =0; i<m; i++){
		A[i] = A[i] + 4.64*static_cast<float>(i)/m;
			}


// the inlet boundary condition
vector<float> rho1(nsp,1e-3);
float u1   = 3500;
float p1   = 2620;
float Tv1  =  500;

rho1[0] = 2e-4;
rho1[1] = 5e-4;
rho1[2] = 3e-4;

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


	for (int k = 0; k<t-1; k++){
		// inlet condition
		for (int j=0; j<nsp; j++){
		q(j,0,k) = rho1[j];	
		}
		q(U,0,k) =   u1;	
		q(P,0,k) =   p1;                              

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
		}
	}
	write_cl( q, m, x, t, nsp, wsp, "../out/cl.dat");


  	return 0;
}


