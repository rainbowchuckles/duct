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

extern "C" {
    void dgetrf_(int*, int*, double*, int*, int*, int*);
    void dgetrs_(char*, int*, int*, double*, int*,
                 int*, double*, int*, int*);
}

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

    double q[NSPMAX],
    double f[NSPMAX],
    double J[NSPMAX][NSPMAX]
);

extern "C" void s2c_c(
    double qp[NSPMAX],
    double qc[NSPMAX],
    double j[NSPMAX][NSPMAX],
    int *nsp,
    double wsp[NSPMAX],
    double asp[NSPMAX][NTRMAX][NACOEF],
    double hfsp[NSPMAX],
    double rsp[NSPMAX][4]
);

extern "C" void o2c_c(
    double qp[NSPMAX],
    double qc[NSPMAX],
    double j[NSPMAX][NSPMAX],
    int *nsp,
    double wsp[NSPMAX],
    double asp[NSPMAX][NTRMAX][NACOEF],
    double hfsp[NSPMAX],
    double rsp[NSPMAX][4]
);

extern "C" void snd_c(
    int *nsp,
    double wsp[NSPMAX],
    double asp[NSPMAX][NTRMAX][NACOEF],
    double hfsp[NSPMAX],
    double rsp[NSPMAX][4],
    double qp[NSPMAX],
    double *a
);


// flux function

void flux(
    double SL[NSPMAX],
    double SR[NSPMAX],
    double F[NSPMAX],
    //double A[NSPMAX][NSPMAX],
    int *nsp,
    double wsp[NSPMAX],
    double asp[NSPMAX][NTRMAX][NACOEF],
    double hfsp[NSPMAX],
    double rsp[NSPMAX][4]
	) {

// left and right fluxes
double A[NSPMAX][NSPMAX];
double FL[NSPMAX] = {0.0};
double FR[NSPMAX] = {0.0};
double UL[NSPMAX] = {0.0};
double UR[NSPMAX] = {0.0};

// nv is the number of variables
// nsp species + Tv + u + p 
int nv  = *nsp + 3;
int Vs  = *nsp;
int Us  = *nsp+1;
int Ps  = *nsp+2;

int I = *nsp;
int E = *nsp+1;
int X = *nsp+3;

// speed of sound
double a;
double aL;
double aR;

// 1. Evaluate U(S) = [rhoi, ..., rhos, rhoev,rhoE,rhou] -> vector of conserved variables

s2c_c(SL, UL, A, nsp, wsp, asp, hfsp, rsp);
s2c_c(SR, UR, A, nsp, wsp, asp, hfsp, rsp);

// 2. Evaluate F(S) = [rhoiu, ..., rhosu, rhoevu,rhoHu,rhou^2 + p] -> flux of conserved variables

FL[E] = UL[E] + SL[Ps];
FR[E] = UR[E] + SR[Ps];

for (int j=0; j<X+1; j++){
	FL[j] = UL[j]*SL[Us];
	FR[j] = UR[j]*SR[Us];
}

FL[X] += SL[Ps];
FR[X] += SR[Ps];

// 6. For each streamwise cell i evaluate the numerical fluxes:
// F_{i+1/2} = 0.5 * ( F(S_L) + F(S_R) )
//             - 0.5 * alpha_{i+1/2} * ( U(S_R) - U(S_L) )
//
// where:
//   S_L  = state to the left of the interface
//   S_R  = state to the right of the interface
//   F(.) = physical flux function
//   U(.) = conserved variable vector
//   alpha_{i+1/2} = max wave speed at the interface

snd_c(nsp,wsp,asp,hfsp,rsp,SL,&aL);
snd_c(nsp,wsp,asp,hfsp,rsp,SR,&aR);

a = max(abs(SL[Us])+aL,abs(SR[Us])+aR);

for (int i=0; i<X+1; i++){
	F[i]  =  0.5*(FL[i]+FR[i]);
	F[i] += -0.5*a*(UR[i]-UL[i]);
}
}


// flux function

void flux1(
    double SL[NSPMAX],
    double F[NSPMAX],
    //double A[NSPMAX][NSPMAX],
    int *nsp,
    double wsp[NSPMAX],
    double asp[NSPMAX][NTRMAX][NACOEF],
    double hfsp[NSPMAX],
    double rsp[NSPMAX][4]
	) {

// left and right fluxes
double A[NSPMAX][NSPMAX];
double FL[NSPMAX] = {0.0};
double FR[NSPMAX] = {0.0};
double UL[NSPMAX] = {0.0};
double UR[NSPMAX] = {0.0};

// nv is the number of variables
// nsp species + Tv + u + p 
int nv  = *nsp + 3;
int Vs  = *nsp;
int Us  = *nsp+1;
int Ps  = *nsp+2;

int I = *nsp;
int E = *nsp+1;
int X = *nsp+3;

// speed of sound
double a;
double aL;
double aR;

// 1. Evaluate U(S) = [rhoi, ..., rhos, rhoev,rhoE,rhou] -> vector of conserved variables

s2c_c(SL, UL, A, nsp, wsp, asp, hfsp, rsp);

// 2. Evaluate F(S) = [rhoiu, ..., rhosu, rhoevu,rhoHu,rhou^2 + p] -> flux of conserved variables

FL[E] = UL[E] + SL[Ps];

for (int j=0; j<X+1; j++){
	FL[j] = UL[j]*SL[Us];
}

FL[X] += SL[Ps];

for (int j=0; j<X+1; j++){
	F[j] = FL[j];
}
}

// robust linear interpolation
void linterp(const double* x1,
             const double* y1,
             int l1,
             const double* x2,
             double* y2,
             int l2)
{
    if (!x1 || !y1 || !x2 || !y2) return;
    if (l1 < 2 || l2 <= 0) return;

    int i = 0;  // interval index for x1

    for (int j = 0; j < l2; ++j)
    {
        double xx = x2[j];

        // ---- BELOW RANGE ----
        if (xx <= x1[0])
        {
            y2[j] = y1[0];
            continue;
        }

        // ---- ABOVE RANGE ----
        if (xx >= x1[l1 - 1])
        {
            y2[j] = y1[l1 - 1];
            continue;
        }

        // ---- FIND INTERVAL ----
        while (i + 1 < l1 && xx > x1[i + 1])
        {
            ++i;
        }

        // Linear interpolation inside interval [i, i+1]
        double xL = x1[i];
        double xR = x1[i + 1];
        double yL = y1[i];
        double yR = y1[i + 1];

        double t = (xx - xL) / (xR - xL);
        y2[j] = yL + t * (yR - yL);
    }
}

// extract column

std::vector<double> getColumn(const std::vector<std::vector<double>>& data,
                              size_t columnIndex)
{
    std::vector<double> column;

    for (const auto& row : data)
    {
        if (columnIndex < row.size())
            column.push_back(row[columnIndex]);
    }

    return column;
}
