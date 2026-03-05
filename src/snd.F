c ... wrapper to get speed of sound
      recursive subroutine snd_c(nsp_c, wsp_c, asp_c, hfsp_c, 
     2                          rsp_c, qp_c, a_c)
     3                bind(C, name="snd_c")
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

      real*8    qw(NSPMAX)
      real*8    awp(NSPMAX,NSPMAX)
      real*8    da2(NSPMAX)

c ... executable statements

      call s2o(qp_c, qw, awp, nsp_c, wsp_c)
      call a2_gnoffo( nsp_c, rsp_c,asp_c,hfsp_c,qw, a_c,da2 )
      a_c = sqrt(a_c)

      return

      end
