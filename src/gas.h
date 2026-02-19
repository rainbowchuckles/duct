#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;

float calp(
    double S[NSPMAX],
    int nsp,
    double wsp[NSPMAX] 
) {
	float p = 0; // working variable
		     //
	for (int j = 1; j<nsp; j++){p += S[j]*8.314*S[nsp+1]/wsp[j];} 

 	p += S[0]*8.314*S[nsp]/wsp[0]; 	

	return p;
}

float calt(
    double S[NSPMAX],
    int nsp,
    double wsp[NSPMAX] 
) {
	float sum = 0; // working variable
	for (int j = 1; j<nsp; j++){sum += S[j]/wsp[j];} 
	// calculate tt using boyle's law
	float tt  = S[nsp+2];
	tt -= S[0]*8.314*S[nsp]/wsp[0];
        tt /= 8.314*sum;	
	
	return tt;
}
