extern "C" void tcw_c(
    const char* str,
//    const int len,
    int *nsp,

    double nelsp[NSPMAX],

    // Fortran: ielsp(NELMAX, NSPMAX)
    double ielsp[NSPMAX][NELMAX],

    // Fortran: melsp(NELMAX, NSPMAX)
    double melsp[NSPMAX][NELMAX],

    double wsp[NSPMAX],

    // Fortran: rsp(4, NSPMAX)
    double rsp[NSPMAX][4],

    // Fortran: asp(NACOEF, NTRMAX, NSPMAX)
    double asp[NSPMAX][NTRMAX][NACOEF],

    double hfsp[NSPMAX],

    // Fortran: mw(3, NSPMAX, NSPMAX)
    double mw[NSPMAX][NSPMAX][3],

    // Fortran: cs(2, NSPMAX, NSPMAX)
    double cs[NSPMAX][NSPMAX][2],

    double *diss,
    double *inz,

    // Fortran: apb(3, NSPMAX)
    double apb[NSPMAX][3],

    int nrn[2],

    // Fortran: nsprn(4, NRNMAX)
    int nsprn[NRNMAX][4],

    // Fortran: isprn(NSPMAX, NRNMAX)
    int isprn[NRNMAX][NSPMAX],

    // Fortran: msprn(NSPMAX, NRNMAX)
    double msprn[NRNMAX][NSPMAX],

    // Fortran: ktbrn(NRNMAX)
    int ktbrn[NRNMAX],

    // Fortran: xtbrn(NSPMAX, NRNMAX)
    double xtbrn[NRNMAX][NSPMAX],

    // Fortran: arr(64, NSPMAX)
    double arr[NRNMAX][64]

);

extern "C" void src_c(
    int *nsp,

    double nelsp[NSPMAX],

    double ielsp[NSPMAX][NELMAX],
    double melsp[NSPMAX][NELMAX],

    double wsp[NSPMAX],

    double rsp[NSPMAX][4],
    double asp[NSPMAX][NTRMAX][NACOEF],

    double hfsp[NSPMAX],

    double mw[NSPMAX][NSPMAX][3],
    double cs[NSPMAX][NSPMAX][2],

    double *diss,
    double *inz,

    double apb[NSPMAX][3],

    int nrn[2],
    int nsprn[NRNMAX][4],
    int isprn[NRNMAX][NSPMAX],
    double msprn[NRNMAX][NSPMAX],

    int ktbrn[NRNMAX],
    double xtbrn[NRNMAX][NSPMAX],

    double arr[NRNMAX][64],

    double *q,
    double *f
);
