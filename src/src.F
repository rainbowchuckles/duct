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

      integer           i,k,flg

      real*8            dt,t_eq,rd

      real*8            ht
         
      real*8            tau
      dimension         tau(NSPMAX,NSPMAX)
      
      real*8            df0
      dimension         df0(NSPMAX,NSPMAX)

      real*8            df1
      dimension         df1(NSPMAX,NSPMAX)

      real*8            dg0
      dimension         dg0(NSPMAX,NSPMAX)

      real*8            dg1
      dimension         dg1(NSPMAX,NSPMAX)

      real*8            dr
      dimension         dr(NSPMAX)
      
      real*8            g0
      dimension         g0(NSPMAX)

c ... executable statements
      frlx = 1.d0
      call src( nsp_c,nelsp_c,ielsp_c,melsp_c,wsp_c,rsp_c,
     1          asp_c,hfsp_c,mw_c,cs_c,diss_c,inz_c,apb_c,
     1          nrn_c,nsprn_c,isprn_c,msprn_c,ktbrn_c,xtbrn_c, 
     1          arr_c,t_eq,qp_c,qp_c,rd,f_c,df0,df1,g0,dg0,
     1          dg1,dr,i,k,
     1          dt,qp_c,flg, .true., 0.35d0,frlx,tau)


      return 

      end

