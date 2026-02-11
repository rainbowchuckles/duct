c ... this subroutine is a wrapper for the OCEAN chmf and thrm files  
c     that initialise the thermochemistry 
      
      subroutine gam_c(wsp_c,asp_c,hfsp_c,
     3                 rsp_c,nsp_c,qp_c,c_c) bind(C, name="gam_c")

      use iso_c_binding, only: c_int, c_double, c_char, c_null_char

      implicit none

c ... include

#     include <typd.h>
#     include <size.h>
#     include <prop.h>
#     include <vars.h>
#     include <park/sizes.h>
#     include <parkvars.h>
#     include <def.h>
#     include <cvars.h>

c ... name of the thermo namelist file

c ... arguments

      character*NFL       trnf
      character*NFL       cond
      character*NFL       chmf
      character*NFL       mw_file
      character*NFL       cs_file
      character*NFL       mst_file
      character*NFL       diss_file
      character*NFL       ion_file
      character*NFL       apb_file
      character*NFL       thrf
      character*NFL       colpth
      character*NFL       rcnf  

c ... file name list

      namelist /chem/  chmf, rcnf,
     1                 mw_file, cs_file, mst_file, diss_file, ion_file,
     2                 apb_file, thrf, colpth!, non_pref, c1

c ... local symbols

      integer           i, l, iu, ios, n,nv

      logical           kill, res
      real*8            c
      real*8            qs 
      dimension         qs(M)
      real*8            q  
      dimension         q(M)
      real*8            dc  
      dimension         dc(M)

c ... executable statements

c ... read reactions file(s)
      nv_c = nsp + 3 

      call snd(nsp_c,nv_c,qp_c,qp_c,c_c,dc,wsp_c,asp_c,hfsp_c,rsp_c)

      return

      end

