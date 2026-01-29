      REAL_                  pi
      parameter            ( pi= 3.14159265358979323844d0 )


         REAL_           runi,patm,kbz,debye,nav,
     1                   trace, logtrace,eps0,ev,
     1                   ec,sbz
         integer         neps
         parameter     ( kbz=1.380649d-23,
     1                   nav=6.02214076d23, 
     2                   runi=nav*kbz,
     3                   patm=101325d0,
     5                   debye=(1.d-21)/299792458d0,
     6                   trace=1.d-20 ,
     7                   logtrace=1.0d-16,
     8                   neps=5,
     9                   eps0=8.8541878128d-12)
         parameter      (ev=1.60217662d-19,
     1                   ec=4.803 204 251d-10 ,
     2                   sbz= 5.670 374 419d-8 )

#ifdef _C_
         const REAL_     kbz=1.380649e-23;
         const REAL_     nav=6.02214076e23;
         const REAL_     runi=nav*kbz;
         const REAL_     patm=101325e0;
         const REAL_     debye=(1.e-21)/299792458e0;
         const REAL_     trace=1.e-20;
         const REAL_     logtrace=1.e-16;
         const REAL_     eps0=8.8541878128e-12
         const REAL_     ev=1.60217662e-19
         const REAL_     ec=4.803 204 251e-10
         const REAL_     sbz= 5.670 374 419e-8
         const int       neps=5;
#endif
