
      extern "C" void init_( int *nel,char *cel, double *wel,
                             int *ise, int *nsp, char *csp,
                             int *nes, int *ies, int *mes,
                             double *wsp, int *ntr, double *mst, 
                             double *amw, double *apb,  
                             int *nrn, int *nsprn, int *isprn, 
                             int *psprn, double *msprn, 
                             int *ntbrn, int *ktbrn, int *itbrn, 
                             double *xtbrn, int *keqrn,
                             double *arn, double *ern, double *trn,
                             int *mspel,
                             int *jspel, int *refsp, int *crsp, 
                             int *nrnsp,int *irnsp,int *prnsp,
                             double *dmrn,int *spar,int *ngs,
                             double *asp,double *rsp,double *hfsp,
                             double *omg_tab,int *omg_ndat,
                             double *coll_inf,double *sc_att,
                             double *sc_rep,int *sc_att_n,
                             int *sc_rep_n,int *chrge_flg);

      extern "C" void thrd_( int *nel, char *cel, int *ise, int *nsp,
                             char *csp, int *nes, int *ies, int *mes,
                             double *rsp, double *asp );

      extern "C" void rcrd_( int *nsp, char *csp, double *wsp, int *nes,
                             int *ies, int *mes, int *nrn, 
                             double *msr, int *nsr, int *isr, int *itr, 
                             int *ntr, double *arr );

      extern "C" void trnrd_( int *nsp,char *csp, const char *name, int *m,
                              double *dij, const char *lbl );

      extern "C" void shk_( int *nsp, double *rsp, double *asp,  double *qi,
                            double *qe, double *dp, int *flag );
      
      extern "C" void  slv_( int *nsp, int *nte, int *nel, char *csp,
                             int *nes, int *ies, int *mes, int *mspel, 
                             int *ise, double *wsp, int *refsp,
                             double *rsp, double *asp, double *hfsp,
                             double *amw, double *apb,
                             int *nrn, int *nsr, int *isr, double *msr,  
                             int *itr, double *mst, double *arr, 
                             int *rig_trans, double *omg_tab, int *omg_ndat,
                             double *coll_inf, double *sc_att, double *sc_rep,
                             int *sc_att_n, int *sc_rep_n, int *chrge_flg,
                             int *n, double *y, double *qi, double *qp, double *dp, 
                             double *b, double *del, double *q, 
                             double *rlx,int *nit, int *flag, double *eps );

      extern "C" void out_( int *nsp, char *csp, double *wsp, double *rsp,
                            double *asp, double *hfsp, int *mm, int *n, 
                            double *b, double *dp, double *y, 
                            double *del, double *qi, double *q, int *l );

      extern "C" void ref_( int *nsp, int *n, int *irf, double *y,
                            double *q, double *c, double *w  );

//      extern "C" void restart_( int *nsp, int *n, double *b, double *dp, double *y, double *del, double *qi, double *q, char *name );
//      extern "C" void save_( int *nsp, int *n, double *b, double *dp, double *y, double *del, double *qi, double *q, char *name );
      extern "C" void restart_( int *l, int *mm, int *n, double *b, double *dp, double *y, double *del, double *qi, double *q, char *name );
      extern "C" void save_( int *l, int*mm, int *n, double *b, double *dp, double *y, double *del, double *qi, double *q, char *name );
