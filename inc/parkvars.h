
      integer               nel
      character(NEL)        cel
      dimension             cel(NELMAX)

      real*8                wel
      dimension             wel(NELMAX)


      integer               nsp, nte
      character(NSP)        csp
      dimension             csp(NSPMAX)


c ... Number of reactions
      integer               ntr
      integer               nrn
      dimension             nrn(2)



c ... Arrhenius constants

      real*8                arr
      dimension             arr(64,NRNMAX)




c ... Molecular weights

      real*8                wsp
      dimension             wsp(NSPMAX)

c ... Specific heat coefficients
c ... temperature ranges

      real*8           rsp
      dimension       rsp(4,NSPMAX)

      real*8                asp
      dimension             asp(NACOEF,NTRMAX,NSPMAX)

      real*8          hfsp
      dimension       hfsp(NSPMAX)

c ... Coefficients for Millikan and White temperature relaxation

      real*8                mw,cs,diss,inz
      dimension             mw(3,NSPMAX,NSPMAX),cs(2,NSPMAX,NSPMAX),
     1                      diss(NSPMAX),inz(NSPMAX)

c ... Polynomial fit for electron neutral energy exchange cross section

      real*8                apb
      dimension             apb(3,NSPMAX)


      integer         refsp,mspel
      dimension       refsp(NSPMAX),mspel(NELMAX,NSPMAX)

c ... correspondence between elements and atomic species

      integer         ispel,        jspel
      dimension       ispel(NELMAX),jspel(NSPMAX)

      integer         nelsp,melsp,ielsp
      dimension       nelsp(NSPMAX),ielsp(NELMAX,NSPMAX),
     1                melsp(NELMAX,NSPMAX)

c ... Rearranging the list of species so electrons are first
c ... and condensed species are last

      integer         spar, ngs
      dimension       spar(NSPMAX)


c ... charge count

      integer         crsp
      dimension       crsp(NSPMAX)

c ... Transport coefficients

      real*8                dij,       bij,
     1                      o11,       o22
      dimension             dij(NOCOEF,NTRMAX,NPRMAX), 
     1                      bij(NOCOEF,NTRMAX,NPRMAX),
     2                      o11(NOCOEF,NTRMAX,NPRMAX),
     3                      o22(NOCOEF,NTRMAX,NPRMAX)

c ... Collision integral and transport data
      real*8          omg_tab,coll_inf(NSPMAX*(NSPMAX+1)/2)
      integer         omg_ndat(NSPMAX*(NSPMAX+1)/2)
      real*8          omg_dat2
      dimension       omg_dat2(11,1024,64*(64+1)/2)
      dimension       omg_tab(11,NLINE,NSPMAX*(NSPMAX+1)/2)
      real*8          sc_att,sc_rep
      dimension       sc_att(11,NLINE),sc_rep(11,NLINE)
      integer         sc_att_n,sc_rep_n,chrge_flg
      real*8          omg_dat
      dimension       omg_dat(10,NSPMAX*(NSPMAX+1)/2)

c ... reactions


      integer         nsprn,        isprn,
     1                psprn
      dimension       nsprn(4,NRNMAX),isprn(NSPMAX,NRNMAX) ,
     1                psprn(NRNMAX)

      integer         ntbrn,ktbrn
      dimension       ntbrn(NRNMAX),ktbrn(NRNMAX)

      integer         keqrn
      dimension       keqrn(NRNMAX)

      integer         itbrn
      dimension       itbrn(NSPMAX,NRNMAX)


      real*8           xtbrn
      dimension       xtbrn(NSPMAX,NRNMAX)

      real*8           arn,ern,trn
      dimension       arn(NRNMAX),ern(NRNMAX),trn(NRNMAX)

      integer         nrnsp,        irnsp,
     1                prnsp
      dimension       nrnsp(NSPMAX),irnsp(NRNMAX,NSPMAX) ,
     1                prnsp(NSPMAX)

      real*8           msprn,                dmrn
      dimension       msprn(NSPMAX,NRNMAX), dmrn(NRNMAX)
