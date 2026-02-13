#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;


float calt(
    const FlowField& q,
    int i,
    int k,
    int nsp,
    double wsp[NSPMAX] 
) {
	float sum = 0; // working variable
	for (int j = 1; j<nsp; j++){sum += q(j,i,k)/wsp[j];} 
	// calculate tt using boyle's law
	float tt  = q(nsp+2,i,k);
	tt -= q(0,i,k)*8.314*q(nsp,i,k)/wsp[0];
        tt /= 8.314*sum;	
	
	return tt;
}

float rmix(
    const FlowField& q,
    int i,
    int k,
    int nsp,
    double wsp[NSPMAX] 
) {
	float sum = 0; // working variable
        float rmix = 0;

	for (int j=0; j<nsp; j++){sum += q(j,i,k)/wsp[j];}
        rmix = 8.314/sum;		
	return rmix;
}
