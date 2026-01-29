#  ifndef _SLINE_
#  define _SLINE_

#  include <cstdlib>
#  include <cstdio>
#  include <cassert>
#  include <cstring>

#  include <park/sizes.h>
#  include <size.h>
#  include <fort.h>
#  include <park.h>

   struct ctrl_t
  {
      double       rlx0;
      double       drlx;
      int          n;
      int          flag;
      int          geo;

      void safe()
     {
         rlx0= 0.1;
         drlx= 1.2;
         n= 20;
         flag= 1;
     };
      void fast()
     {
         rlx0= 1.0;
         drlx= 1.2;
         n= 10;
         flag= 0;
     }
  };

   struct sline_t
  {
      protected:

         double                dp;
         double                qe[NSPMAX];

// elemental species

         int                   nte;

         int                   nel;
         char                  cel[NEL*NELMAX];

         double                wel[NELMAX];

         int                   ise[NELMAX];

         int                   nsp;
         char                  csp[NSP*NSPMAX];


//... chemical species compositions

         int                   nes[NSPMAX],ies[NELMAX*NSPMAX],
                               mes[NELMAX*NSPMAX];

//... Number of reactions

         int                   nrn[2];

//... Number of species per reaction

         int                   nsr[4*NRNMAX];

//... Stoichiometric coefficients

         double                msr[64*NRNMAX];

//... Species per reaction

         int                   isr[64*NRNMAX];

//... Arrhenius constants

         double                arr[64*NRNMAX];

//... Number of third-body types

         int                   ntr;
         int                   itr[NRNMAX];

//... Third-body coefficients

         double                mst[NSPMAX*NSPMAX];

//... Molecular weights

         double                wsp[NSPMAX];

//... Specific heat coefficients

         double                asp[NACOEF*NTRMAX*NSPMAX];

//... Coefficients for Millikan and White temperature relaxation

         double                amw[4*NSPMAX];

//... Polynomial fit for electron neutral energy exchange cross section

         double                apb[3*NSPMAX];

//... fits for collision integrals and their ratios

         double                dij[NOCOEF*NTRMAX*NPRMAX], 
                               bij[NOCOEF*NTRMAX*NPRMAX],
                               o11[NOCOEF*NTRMAX*NPRMAX],
                               o22[NOCOEF*NTRMAX*NPRMAX];

//... for reading chem and thermo data
         int                   nsprn[NRNMAX],
                               isprn[NSPMAX*NRNMAX],
                               psprn[NRNMAX];

         int                   ntbrn[NRNMAX],
                               ktbrn[NRNMAX];
 
         int                   keqrn[NRNMAX];

         int                   itbrn[NSPMAX*NRNMAX];
   
   
         double                xtbrn[NSPMAX*NRNMAX];
   
         int                   nrnsp[NSPMAX],
                               irnsp[NRNMAX*NSPMAX],
                               prnsp[NSPMAX];
   
         double                msprn[NSPMAX*NRNMAX], 
                               dmrn[NRNMAX];

         double                arn[NRNMAX],
                               ern[NRNMAX],
                               trn[NRNMAX];
   
         int                   ngs;

         int                   spar[NSPMAX];
// ... temperature ranges
   
         double                rsp[4*NSPMAX];
  
   
         double                hfsp[NSPMAX];
   
// ... correspondence between elements and atomic species
   
         int                  ispel[NELMAX],
                              jspel[NSPMAX];
   
         int                  refsp[NSPMAX],
                              mspel[NELMAX*NSPMAX];
   
         
// ...    charge count
   
         int                  crsp[NSPMAX];
   
              
// ...    Collision integral and transport data
         double               coll_inf[(NSPMAX*(NSPMAX+1)/2)];

         int                  omg_ndat[(NSPMAX*(NSPMAX+1)/2)];

         double               omg_tab[11*NLINE*(NSPMAX*(NSPMAX+1)/2)];

         double               sc_att[11*NLINE],
                              sc_rep[11*NLINE];

         int                  sc_att_n,
                              sc_rep_n,
                              chrge_flg;

         double               omg_dat[10*(NSPMAX*(NSPMAX+1)/2)];

         int                  rig_trans;

     public:
        void init();

        void    atm( int &m, double *q);

        void   solve( double *qi, double *qp, double b, int *n, double *y, double *q, double &del, ctrl_t c, double eps );

        void    save( double *qi, double b, int *n, double *y, double *q, double &del, char *name );

        void restart( double *qi, double &b, int *n, double *y, double *q, double &del, char *name );

       void    guess( double *qi, int n, double *q );
       void    guess( double *qi, int n0, int n1, double *q );
  
       void   refine( double *qi, int &n, double *y, double *q );
  };

#endif
