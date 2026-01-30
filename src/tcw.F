c ... this subroutine is a wrapper for the OCEAN chmf and thrm files  
c     that initialise the thermochemistry 
      
      subroutine tcw_c(nsp_c,nelsp_c,ielsp_c,melsp_c,wsp_c,rsp_c,
     1                 asp_c,hfsp_c,mw_c,cs_c,diss_c,inz_c,apb_c,
     2                 nrn_c,nsprn_c,isprn_c,msprn_c,ktbrn_c,xtbrn_c,
     3                 arr_c) bind(C, name="tcw_c")

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

      character*NFL       thrf,colpth
      character*NFL       chmf
      character*NFL       rcnf
      character*NFL       trnf
      character*NFL       cond
      character*NFL       mw_file
      character*NFL       cs_file
      character*NFL       diss_file
      character*NFL       ion_file
      character*NFL       apb_file
      character*NFL       mst_file

c ... local symbols

c ... file paths
      chmf     = '../../slime/in/thermo/park93-6sp/chemkin_o2_n2.dat' 
      rcnf     = '../../slime/in/thermo/park93-6sp/'//
     &           'reactions_park93_o2_n2.dat' 
      mw_file  = '../../slime/in/thermo/park93-6sp/'//
     &           'mw_millikan_o2_n2.dat' 
      cs_file  = '../../slime/in/thermo/park93-6sp/cs_o2_n2.dat' 
      mst_file = '../../slime/in/thermo/park93-6sp/tbr_park93_o2_n2.dat' 
      diss_file= '../../slime/in/thermo/park93-6sp/diss_o2_n2.dat' 
      ion_file = '../../slime/in/thermo/park93-6sp/'//
     &           'ion_gnoffo89_o2_n2.dat' 
      apb_file = '../../slime/in/thermo/park93-6sp/e_neu_cs_o2_n2.dat'
      thrf     = '../../slime/in/thermo/thermo.inp' 
      colpth   = '../../slime/in/thermo/bin/omg_beta/'

c ... executable statements

      call  chm( chmf, nel, cel, wel,
     1           nsp, csp,
     1           nrn, nsprn, isprn, psprn, msprn,
     1           ntbrn,  ktbrn,  itbrn,  xtbrn,  keqrn,
     1           arn, ern, trn  )

      call thrm( thrf, nel, cel, wel,
     1           nsp, csp, nelsp, melsp, mspel,ielsp,
     1           ispel, jspel, refsp, crsp, wsp,
     1           nrn, nsprn,isprn,psprn,nrnsp,irnsp,prnsp,
     1           msprn,dmrn,spar,ngs,asp,rsp,hfsp, .false.)

c ... initialize non-eq model

      call init( nsp,mw_file,cs_file,diss_file,ion_file,apb_file,
     1       mst_file,ntr,xtbrn, mw,cs,diss,inz,apb )

c ... read thermo.inp

      call thrd( nel,cel,ispel, nsp,csp,nelsp,ielsp,melsp,rsp,asp)

c ... read reactions file(s)

      print *, ' '
      print *, 'Reaction scheme: '
      print *, ' '

      call rcrd(  nsp,csp,wsp, nelsp,ielsp,melsp,
     1                 nrn, msprn,nsprn,isprn,ktbrn,ntr,arr ,rcnf)

      print *, ' '
      print *, '---------------------'
      print *, ' '

      nsp_c = nsp
      nelsp_c = nelsp
      ielsp_c = ielsp
      melsp_c = melsp
      wsp_c = wsp  
      rsp_c = rsp  
      asp_c = asp  
      hfsp_c = hfsp  
      mw_c = mw  
      cs_c = cs  
      diss_c = diss  
      inz_c = inz  
      apb_c = apb  
      nrn_c = nrn  
      nsprn_c = nsprn  
      isprn_c = isprn  
      msprn_c = msprn  
      ktbrn_c = ktbrn  
      xtbrn_c = xtbrn  
      arr_c = arr  

      return

      end

