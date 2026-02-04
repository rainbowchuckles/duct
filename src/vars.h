// parameters for the non-equilibrium model
// see OCEAN for definitions 
double nelsp[NSPMAX];
double ielsp[NSPMAX][NELMAX];
double melsp[NSPMAX][NELMAX];
double wsp[NSPMAX]; // Fix
double rsp[NSPMAX][4];
double asp[NSPMAX][NTRMAX][NACOEF];
double hfsp[NSPMAX];
double mw[NSPMAX][NSPMAX][3];
double cs[NSPMAX][NSPMAX][2];
double diss[NSPMAX];
double inz[NSPMAX];
double apb[NSPMAX][3];
int    nrn[2];
int    nsprn[NRNMAX][4];
int    isprn[NRNMAX][NSPMAX];
int    msprn[NRNMAX][NSPMAX];
int    ktbrn[NRNMAX];
int    xtbrn[NRNMAX][NSPMAX];
double arr[NSPMAX][64];

// number of species
int nsp;

// qp is the vector of primitive variables
// q(0:nsp-1) = partial densities of species, kg/m3
// q(V)       = vibrational-electronic temperature, K
// q(U)       = velocity in the lab frame, m/s
// q(P)       = pressure, Pa
double qp[NSPMAX];

// f is the vector of non-equilibrium source terms from OCEAN/LASTA
// they are returned as /s so need to be multiplied by dt before 
// accumulating with the flow vector
double f[NSPMAX];

// gas
// the ratio of specific heats
float gam;

// time       
// t is the number of time steps
int t;

// dt is the size of the time step
float dt;

// grid
// m is the number of grid points
int m;

// l is the length of the domain
float l;

// the inlet boundary condition
float u1;
float p1;
float Tv1;

// solution controls
float conv;

// thermochemistry input files
string chmf;
string rcnf; 
string mw_file;
string cs_file; 
string mst_file; 
string diss_file; 
string ion_file; 
string apb_file; 
string thrf; 
string colpth;
