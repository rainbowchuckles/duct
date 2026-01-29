cc ... 
c      integer        nsp

c ... thermochemical model
c
c      REAL_          rsp
c      dimension       rsp(4,64)

c ... standard formation properties 

c      REAL_          asp
c      dimension       asp(10,3,64)

c      REAL_          hfsp
c      dimension       hfsp(64)
c ... jspdec

      !integer         nsp,ngs
c      integer          ngs

c      character*16    csp
c      dimension       csp(64)

c ... Elemental compositions of molecular species

c     nelsp           number of elements present in each species
c     ielsp           identiy of the elements present in each species
c                     in the order read from the chemkin.dat file
c     melsp           amount of the each element presen in the species

c      integer         nelsp,melsp,ielsp,refsp,mspel
c      dimension       nelsp(64),ielsp(8,64),
c     1                melsp(8,64),
c     1                refsp(64),mspel(8,64)

c ... correspondence between elements and atomic species

c      integer         ispel,        jspel
c      dimension       ispel(8),jspel(64)

c ... Rearranging the list of species so electrons are first
c ... and condensed species are last

c      integer         spar
c      dimension       spar(64)
c ... molecular masses
c
c      REAL_           wsp
c      dimension       wsp(64)

c ... charge count
c
c      integer         crsp
c      dimension       crsp(64)

c ... jrn dec
c      integer         nrn

c      integer         nsprn,        isprn,
c     1                psprn
c      dimension       nsprn(128),isprn(64,128) ,
c     1                psprn(128)

c      integer         ntbrn,ktbrn
c      dimension       ntbrn(128),ktbrn(128)

c      integer         keqrn
c      dimension       keqrn(128)

c      integer         itbrn
c      dimension       itbrn(64,128)

c      REAL_           xtbrn
c      dimension       xtbrn(64,128)

c      REAL_           arn,ern,trn
c      dimension       arn(128),ern(128),trn(128)

c      integer         nrnsp,        irnsp,
c     1                prnsp
c      dimension       nrnsp(64),irnsp(128,64) ,
c     1                prnsp(64)
c
c      REAL_           msprn,                dmrn
c      dimension       msprn(64,128), dmrn(128)

c ... el dec
c      integer         nel
c
c      character*2   cel
c      dimension       cel(8)
c
c      REAL_           wel
c      dimension       wel(8)

c ... pack everything from HERE
c ... elements


c ... Elemental compositions of molecular species

c     nelsp           number of elements present in each species
c     ielsp           identiy of the elements present in each species
c                     in the order read from the chemkin.dat file
c     melsp           amount of the each element presen in the species
c     igssp           flag for condensed phases (unused now but
c                     needs to be included when computing partial
c                     pressures)


c ... to HERE in a common and create an equivalent C structure
c     We keep the best of two worlds: the C structure is nicely incapsulated in the FROSST structure,
c     but at the same time is fully accessible from the Fortran functions 
c      REAL_          omg_dat
c     REAL_          omg_dat2
c      REAL_          coll_inf(64*(64+1)/2)!,eta,k_mix
      REAL_          k_e,k_h,k_r,k_int,k_tot,el_cond,dmu
c      integer         omg_ndat(64*(64+1)/2),nval
      dimension       dmu(64+3)
      
c      REAL_          sc_att,sc_rep,v,c,spd
c      dimension       sc_att(11,1024),sc_rep(11,1024)
c      integer         sc_att_n,sc_rep_n,chrge_flg

      REAL_             xin
      dimension          xin(3,64)
      
      REAL_             xout
      dimension          xout(3,64)

      REAL_             exout
      dimension          exout(3,32)

      integer         frz_flg
 




