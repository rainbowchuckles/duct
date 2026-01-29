c ...

c ... Avogadro number

      real*8                nav
      parameter           ( nav=6.022 140 76d23 )
      
c ... Electronvolt value in Joule

      real*8                ev
      parameter           ( ev=1.602 176 634d-19 )

c ... Electron charge in ESU

      real*8                ec
      parameter           ( ec=4.803 204 251d-10 )

c ... Boltzman constant

      real*8                kbz
      parameter           ( kbz=1.380 649 d-23 )

c ... Stefan-Boltzman constant

      real*8                sbz
      parameter           ( sbz= 5.670 374 419 d-8 )

c ... Universal gasl constant
      
      real*8                runi
      parameter           ( runi= kbz*nav )

c ... PI

      real*8                pi
      parameter           ( pi= 3.14159265358979323844d0 )

c ... trace species molality

      real*8                trace
      parameter           ( trace= 1.d-20 )

      real*8                logtrace
      parameter           ( logtrace= 1.d-16 )

      integer               neps
      parameter           ( neps=4 )
