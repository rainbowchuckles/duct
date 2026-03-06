<p align="center">
  <img src="nozzle.png" alt="Project Logo" width="900"/>
</p>


Joseph Steer 6/3/26

This software was developed to model the flow through a supersonic nozzle.

## Governing Equations

A quasi-1D, inviscid, two-temperature variable area duct obeys the conservation equations:

    d/dt (A*U) + d/dx (A*F) = A*Q

where U is the vector of conserved variables, F is the flux vector, Q is the vector of
thermochemical source terms, and A is the area as a function of streamwise distance.

This can be expressed in matrix form as:

         | rho_s |            | rho_s*u         |                | rho_s*u         |   | q_dot_s |
         | ...   |            | ...             |                | ...             |   | ...     |
    d/dt | rho_N | =  -d/dx   | rho_N*u         | - (1/A)(dA/dx) | rho_N*u         | + | q_dot_N |
         | rho*u |            | rho*u*u + p     |                | rho*u*u         |   | 0       |
         |rho*e_V|            | rho*e_V*u + p_e |                | rho*e_V*u + p_e |   | q_dot_V |
         | rho*E |            | rho*H*u         |                | rho*H*u         |   | 0       |

where (1) are the conserved variables, (2) are the fluxes, (3) are the geometric source
terms, and (4) are the thermochemical source terms.

--------------------------------------------------------------------------------

## Numerical Scheme

The system is discretised and solved using a finite-volume scheme with an explicit update:

    S_i^(k+1) = S_i^(k) + dt * inv(J_US)|_i^k
                * [ -(1/dx) * (F_hat_(i+1/2)^k - F_hat_(i-1/2)^k)
                    - (1/A_i^k) * (A_(i+1/2) - A_(i-1/2))/dx * (F_i^k - G_i^k)
                    + Q_i^k ]

where S is the vector of state variables:

    S = | rho_s |
        | ...   |
        | rho_N |
        | u     |
        | T_V   |
        | p     |

--------------------------------------------------------------------------------

The interface flux F_hat_(i+1/2)^k is computed using the Rusanov (local Lax-Friedrichs) scheme:

    F_hat_(i+1/2)^k = 0.5 * (F_i^k + F_(i+1)^k)
                      - 0.5 * alpha_(i+1/2)^k * (U_(i+1)^k - U_i^k)

where alpha is the maximum wave speed at the cell interface:

    alpha_(i+1/2)^k = max( |u_i^k + a_i^k|, |u_(i+1)^k + a_(i+1)^k| )

and u and a are the velocity and speed of sound respectively.


## Quick start

Compile using make, you will need g++, gcc, OCEAN, lapack, lblas, openmpi
Execute using ./a.out air.inp air-6sp-thermo.inp X3-M12.dat
Initial conditions need to be set in main.cpp (will fix this later)

