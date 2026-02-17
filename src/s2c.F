c ... wrapper for lasta routine that converts primitive to conserved variables 

      subroutine s2c_c(qp_c,qc_c, j_c, nsp_c,wsp_c,asp_c,hfsp_c,
     3                 rsp_c) bind(C, name="s2c_c")
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

c ... executable statements

      j = 0.0d0
      qs = 0.0d0
      
      call s2c(qp_c,qc_c,j,nsp_c,wsp_c,asp_c,hfsp_c,rsp_c,qs)

      return 

      end
