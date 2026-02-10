c ... this subroutine is a wrapper for the OCEAN chmf and thrm files  
c     that initialise the thermochemistry 
      
      subroutine tcw_c(inp,nsp_c,nelsp_c,ielsp_c,melsp_c,wsp_c,rsp_c,
     1                 asp_c,hfsp_c,mw_c,cs_c,diss_c,inz_c,apb_c,
     2                 nrn_c,nsprn_c,isprn_c,msprn_c,ktbrn_c,xtbrn_c,
     3                 arr_c) bind(C, name="tcw_c")

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

      character(kind=c_char), dimension(*) :: inp
      character(len=:), allocatable::path

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

      integer             i, l, iu, ios, n

      logical             kill, res

c ... get the file from name list file

c ... find length
      i = 0
      do while (inp(i+1) /= c_null_char)
        i = i + 1
      end do
      
      if (i > 0) then
        allocate(character(len=i) :: path)
        do n = 1, i
          path(n:n) = transfer(inp(n), ' ')   ! convert c_char -> Fortran char
        end do
      else
        path = ''
      end if
      iu = 10
      open(unit=iu, file=path, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,*) 'Input file not found: ', path
         stop
      endif

c ... read

      rewind(iu)
      read(iu, nml=chem,  iostat=ios)
      close(iu)

c ... input validation

      kill = .false.

      inquire( file=trim(chmf), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'chmf file not found: ', chmf
         kill = .true.
      endif

      inquire( file=trim(rcnf), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'rcnf file not found: ', rcnf
         kill = .true.
      endif

      inquire( file=trim(mw_file), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'mw_file file not found: ', mw_file
         kill = .true.
      endif

      inquire( file=trim(cs_file), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'cs_file file not found: ', cs_file
         kill = .true.
      endif

      inquire( file=trim(mst_file), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'mst_file file not found: ', mst_file
         kill = .true.
      endif

      inquire( file=trim(diss_file), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'diss_file file not found: ', diss_file
         kill = .true.
      endif

      inquire( file=trim(ion_file), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'ion_file file not found: ', ion_file
         kill = .true.
      endif

      inquire( file=trim(apb_file), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'apb_file file not found: ', apb_file
         kill = .true.
      endif

      inquire( file=trim(thrf), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'thrf file not found: ', thrf
         kill = .true.
      endif

      inquire( file=trim(colpth), exist=res )
      if (res .eqv. .false.) then
         write(*,*) 'colpth file not found: ', colpth
         kill = .true.
      endif 

c     check if the path provided for colpth is a directory
      l = 0
      do 30 i = len(colpth), 1, -1
      if (colpth(i:i) .ne. ' ') then
            l = i
            go to 40
      end if
30    continue
40    continue

      if (l .gt. 0) then
       if (colpth(l:l) .ne. '/') then
         write(*,*) 'Error: colpth must be a directory'
         stop
       end if
      end if

c ... stop the program if necessary

      if (kill .eqv. .true.) then
         write(*,*) 'Stopping'
         stop
      endif

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

      
      print *, '---------------------'
      print *, ' '
      print *, 'Species: '
      print *, ' '
      do i=1,nsp
      print *, csp(i)
      end do
      print *, ' '
      print *, '---------------------'
      print *, ' '
      print *, 'Reaction scheme: '
      print *, ' '

c ... read reactions file(s)
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

