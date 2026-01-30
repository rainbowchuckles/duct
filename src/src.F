c ... this subroutine is a wrapper for the LASTA/NESS src.F subroutine
c     that calculates non-equilibrium source terms for an arbitrary
c     two-temperature plasma
      
      subroutine src_c(nsp_c,nelsp_c,ielsp_c,melsp_c,wsp_c,rsp_c,
     1                 asp_c,hfsp_c,mw_c,cs_c,diss_c,inz_c,apb_c,
     2                 nrn_c,nsprn_c,isprn_c,msprn_c,ktbrn_c,xtbrn_c,
     3                 arr_c,qp_c,f_c) bind(C, name="src_c")
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

c ... executable statements

      call src( nsp_c,nelsp_c,ielsp_c,melsp_c,wsp_c,rsp_c,
     1          asp_c,hfsp_c,mw_c,cs_c,diss_c,inz_c,apb_c,
     1          nrn_c,nsprn_c,isprn_c,msprn_c,ktbrn_c,xtbrn_c, 
     1          arr_c,0,0,qp_c,0,f_c,0,0,0,0,0,0,0,0,
     1          0,0,0, .false., 0.3, 1.0,0.0)
      return

      end

