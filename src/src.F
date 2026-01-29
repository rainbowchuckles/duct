c ... this subroutine is a wrapper for the LASTA/NESS src.F subroutine
c     that calculates non-equilibrium source terms for an arbitrary
c     two-temperature plasma
      
      subroutine src_c(n) bind(C, name="src_c")

      use iso_c_binding, only: c_int, c_float

      implicit none

      integer(c_int), intent(inout) :: n 

c ... arguments

c ... local symbols

c ... executable statements
        n = n + 3
c      call src( nsp, nelsp,ielsp,melsp, wsp, rsp,asp,hfsp, mw,cs,
c     1          diss,inz,apb,nrn, nsprn,isprn,msprn, ktbrn,xtbrn, 
c     1          arr,t_q, qp1,qp, rd,f,df0,df1,g,dg0,dg1,dr,pi,pk,
c     1          dt,qs,flg, non_pref, c1, frlx,tau)
      
      return

      end

