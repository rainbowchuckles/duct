c ... wrapper for lasta routine that converts primitive to conserved variables 

      subroutine o2c_c(qw_c,qc_c,j, nsp_c,wsp_c,asp_c,hfsp_c,
     3                 rsp_c) bind(C, name="o2c_c")
      use iso_c_binding, only: c_int, c_double

      implicit none

c ... include

#     include <typd.h>
#     include <size.h>
#     include <prop.h>
#     include <vars.h>
#     include <park/sizes.h>
#     include <parkvars.h>
#     include <ocean000.h>
#     include <def.h>
#     include <cvars.h>

c ... arguments

c ... local symbols

      real*8            j(NSPMAX, NSPMAX)

      real*8            qs(NSPMAX)

      real*8            rho

      integer           i

c ... executable statements

      rho = 0.d0

      do 1 i=1,nsp_c
        rho= rho+ qw_c(i)
    1 continue           

c     convert from density to molality

      do 3 i=1,nsp
        qw_c(i) = qw_c(i)/(wsp_c(i)*rho)
    3 continue

      call s2c(qw_c,qc_c,j,nsp_c,wsp_c,asp_c,hfsp_c,rsp_c)

      return 

      end
