c ... pre-defined array sizes

c     M is the maximum number of points in the streamwise direction
c     N is the maximum number of radial positions
c     T is the number of time levels

c 	  L, the number of spatial points
#     define L  601
c 	  O, the number of boundary layer stations (presently unused) 
#     define O  1
c     R, the number of time steps
#     define R  151
c 	  S, size of the pressure profile read from file
#     define S  40000
#     define _N S
c 	  A, the dimension of the training data
#     define A  3 
c 	  B, the number of training points      
#     define B  5000   
c 	  F, the number of candidate points
#     define F  400 

c ... value to be replaced with the macro of same(?) name from OCEAN.
#     define    NSPMX 21 

c ... these are the values that come from NESS/OCEAN

#     define            N 64
#     define            M NSPMAX

#     define            Q 0
#     define            T 1
#     define            P 2
#     define            U 3
#     define            V 4
#     define            W 5

#     define            I 0
#     define            E 1
#     define            C 2
#     define            X 3
#     define            Y 4
#     define            Z 5

#     define            F1 6
#     define            F2 9
