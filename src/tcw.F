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

      chmf     = '../../slime/in/thermo/Park94/chemkin_mars.dat' 
      rcnf     = '../../slime/in/thermo/Park94/reactions_park94.dat' 
      mw_file  = '../../slime/in/thermo/Park94/mw_millikan_park94.dat' 
      cs_file  = '../../slime/in/thermo/Park94/cs_park94.dat' 
      mst_file = '../../slime/in/thermo/Park94/tbr_park94.dat' 
      diss_file= '../../slime/in/thermo/Park94/diss_gnoffo89_park94.dat' 
      ion_file = '../../slime/in/thermo/Park94/ion_gnoffo89_park94.dat' 
      apb_file = '../../slime/in/thermo/Park94/' //
     &           'e_neu_cs_gnoffo89_park94.dat'
      thrf     = '../../slime/in/thermo/thermo.inp' 
      colpth   = '../../slime/in/thermo/bin/omg_dodgy/' 

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

      nelsp_c = nelsp
      ielsp_c = ielsp
      melsp_c = melsp
      wsp_c = wsp  

      return

      end

