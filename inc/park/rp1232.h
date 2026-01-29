c ...
c ... Numerical constants for  NASA RP1232
c
c     R.N. Gupta, J.M. Yos, R.A. Thompson, K.P. Lee, 
c     A review of reaction rates and thermodynamic and transport properties
c     for an 11 species air model for chemical and thermal nonequilibrium
c     calculations to 30000K, NASA RP1232, 1990.
c
c     For details of blending of coefficients across temperature ranges
c     R.A. Thompson, K.P. Lee, R.N. Gupta, Computer codes for the 
c     evaluation of themodynamic properties, transport properties, and
c     equilibrium constants of an 11-species air model.
c     NASA TM 102602, 1990.
c

      real*8                   tref
      parameter              ( tref=298.168d0 )

      real*8                   tr
      dimension                tr(NTRMAX)

c Temperature ranges in NASA RP 1232
c
c     data                     tr( 1) /    300.d0 /,
c    1                         tr( 2) /    800.d0 /,
c    1                         tr( 3) /   1200.d0 /,
c    1                         tr( 4) /   5500.d0 /,
c    1                         tr( 5) /   6500.d0 /,
c    1                         tr( 6) /  14500.d0 /,
c    1                         tr( 7) /  15500.d0 /,
c    1                         tr( 8) /  24500.d0 /,
c    1                         tr( 9) /  25500.d0 /,
c    1                         tr(10) /  30000.d0 /

c Modified temperature ranges to remove oscillations in dcp/dt

      data                     tr( 1) /    300.d0 /,
     1                         tr( 2) /    800.d0 /,
     1                         tr( 3) /   1200.d0 /,
     1                         tr( 4) /   4500.d0 /,
     1                         tr( 5) /   7500.d0 /,
     1                         tr( 6) /  13500.d0 /,
     1                         tr( 7) /  16500.d0 /,
     1                         tr( 8) /  23500.d0 /,
     1                         tr( 9) /  26500.d0 /,
     1                         tr(10) /  30000.d0 /

