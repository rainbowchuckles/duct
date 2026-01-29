c ... Coefficients for Millikan and White temperature relaxation

      real*8                amw
      dimension             amw(4,NSPMAX)

c ... Polynomial fit for electron neutral energy exchange cross section

c      real*8                apb
c      dimension             apb(3,NSPMAX)

      integer               ise
      dimension             ise(NELMAX)  

c ... chemical species compositions

      integer               nes, ies,mes
      dimension             nes(NSPMAX),ies(NELMAX,NSPMAX),
     1                      mes(NELMAX,NSPMAX)

c ... Third-body coefficients

      real*8                mst
      dimension             mst(NSPMAX,NSPMAX)       

c ... Number of third-body types

c      integer               ntr
      integer               itr(NRNMAX)
      integer               itb(NRNMAX)

c ... Collision integral and transport data
c      real*8                omg_tab    

c ... Arrhenius constants

c      real*8                arr
c      dimension             arr(32,NRNMAX)   

c ... Species per reaction

      integer               isr
      dimension             isr(32,NRNMAX)     
      
c ... Stoichiometric coefficients

      real*8                msr
      dimension             msr(32,NRNMAX)     
      
c ... Number of species per reaction

      integer               nsr
      dimension             nsr(4,NRNMAX)   

c ... Third-body coefficients

      real*8                mtb
      dimension             mtb(NSPMAX,NSPMAX)