#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

#include "class.h"
#include "post.h"

using namespace std;

// the main program

int main() {

// gas
// the ratio of specific heats
float gam = 1.4;

// time       
// t is the number of time steps
int t = 15000;
// dt is the size of the time step
float dt = 1.e-7;

// grid
// m is the number of grid points
int m = 400;
// l is the length of the domain
float l = 2.0;
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

// the flow vector
// q(N,M,T)
// N: the primitive variables (rho, u, p)
// M: the streamwise points 
// T: the time steps

FlowField q(3,m,t);

// the initial condition
float rhoi = 1e-3;
float ui   = 3500;
float pi   = 2620;

// the inlet boundary condition
float rho1 = 1e-3;
float u1   = 3500;
float p1   = 2620;

// CFL check
float sgm = ui*(x[2]-x[1])/dt;

	for (int i =0; i<m; i++){
		for (int k= 0; k<t; k++){
			q(0,i,k) = rhoi;
			q(1,i,k) =   ui;
			q(2,i,k) =   pi;
		}
	}


// Euler updates

float dx;
float dA;
float drho;
float du;
float dp;

	for (int k = 0; k<t-1; k++){
		// initial condition
		q(0,0,k) = rho1;	
		q(1,0,k) =   u1;	
		q(2,0,k) =   p1;                              

		for (int i = 1; i<m; i++){
		dx   = x[i] - x[i-1];
		dA   = A[i] - A[i-1];
		drho = q(0,i,k) - q(0,i-1,k);
		du   = q(1,i,k) - q(1,i-1,k);
		dp   = q(2,i,k) - q(2,i-1,k);

		// continuity update
		q(0,i,k+1) = dA/dx;         
		q(0,i,k+1) = q(0,i,k+1)*(q(1,i,k)/A[i]);
		q(0,i,k+1) = q(0,i,k+1) + du/dx;        
		q(0,i,k+1) = q(0,i,k+1)*q(0,i,k);
		q(0,i,k+1) = q(0,i,k+1) + q(1,i,k)*drho/dx; 
		q(0,i,k+1) = dt*q(0,i,k+1);

		q(0,i,k+1) = q(0,i,k) - q(0,i,k+1);

		// momentum update
		q(1,i,k+1) = dp/dx;
		q(1,i,k+1) = q(1,i,k+1)/q(0,i,k);
		q(1,i,k+1) = q(1,i,k+1) + q(1,i,k)*du/dx;
		q(1,i,k+1) = q(1,i,k+1)*dt;              
		
		q(1,i,k+1) = q(1,i,k) - q(1,i,k+1);              

		// energy update
		q(2,i,k+1) = dA/dx;         
		q(2,i,k+1) = q(2,i,k+1)*(q(1,i,k)/A[i]);
		q(2,i,k+1) = q(2,i,k+1) + du/dx;        
		q(2,i,k+1) = q(2,i,k+1)*gam*q(2,i,k);   
		q(2,i,k+1) = q(2,i,k+1) + q(1,i,k)*dp/dx;
		q(2,i,k+1) = q(2,i,k+1)*dt;              
		
		q(2,i,k+1) = q(2,i,k) - q(2,i,k+1);              
		}
	}
	write_flowfield(q, m, t, "../out/flowfield.dat");
	write_outlet(   q, m, t, "../out/outlet.dat");
	write_cl(       q, m, t, "../out/cl.dat");


  	return 0;
}


