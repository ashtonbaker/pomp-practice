glob_snippet <- "
#define ESTAGES 5
#define LSTAGES 7
#define PSTAGES 7
#define ASTAGES 1
#define L_0 250
#define P_0 5
#define A_0 100
"

init_snippet <- "
double *E = &E1;
double *L = &L1;
double *P = &P1;

int k;
double E_tot = 0;
double L_tot = 0;
double P_tot = 0;
for (k = 0; k < ESTAGES; k++) E_tot += E[k];
for (k = 0; k < LSTAGES; k++) L_tot += L[k];
for (k = 0; k < PSTAGES; k++) P_tot += P[k];

double gamma_E = (ESTAGES / tau_E) * exp((-cel * L_tot - cea * A) / ESTAGES);
double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A) / ESTAGES);

double mu_e = (ESTAGES / tau_E) * (1 - exp((-cel * L_tot - cea * A) / ESTAGES));
double mu_l = (LSTAGES / tau_L) * mu_L;
double mu_p = (PSTAGES / tau_P) * (1 - exp((-cpa * A) / ESTAGES));

double L_rate[LSTAGES], P_rate[100] = {0};

for (k = 0; k < LSTAGES; k++) L_rate[k] = pow(gamma_L, k);
for (k = 0; k < PSTAGES; k++) P_rate[k] = pow(gamma_P, k);

for (k = 0; k < ESTAGES; k++) E[k] = 0;
reulermultinom(LSTAGES, L_0, &L_rate[0], 1, &L[0]);
reulermultinom(PSTAGES, P_0, &P_rate[0], 1, &P[0]);
A = 100;"

rproc_snippet <- 
"
double *E = &E1;
double *L = &L1;
double *P = &P1;

int k;
double E_tot = 0;
double L_tot = 0;
double P_tot = 0;
for (k = 0; k < ESTAGES; k++) E_tot += E[k];
for (k = 0; k < LSTAGES; k++) L_tot += L[k];
for (k = 0; k < PSTAGES; k++) P_tot += P[k];

double gamma_E = (ESTAGES / tau_E) * exp((-cel * L_tot - cea * A) / ESTAGES);
double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A) / ESTAGES);

double mu_e = (ESTAGES / tau_E) * (1 - exp((-cel * L_tot - cea * A) / ESTAGES));
double mu_l = (LSTAGES / tau_L) * mu_L;
double mu_p = (PSTAGES / tau_P) * (1 - exp((-cpa * A) / ESTAGES));

double rate[2], etrans[28], ltrans[28], ptrans[28], adeath;

// Calculate who goes where
for (k = 0; k < ESTAGES; k++) {
rate[0] = gamma_E;
rate[1] = mu_e;
reulermultinom(2,E[k],&rate[0],1,&etrans[2*k]);
}

for (k = 0; k < LSTAGES; k++) {
rate[0] = gamma_L;
rate[1] = mu_l;
reulermultinom(2,L[k],&rate[0],1,&ltrans[2*k]);
}

for (k = 0; k < PSTAGES; k++) {
rate[0] = gamma_P;
rate[1] = mu_p;
reulermultinom(2,P[k],&rate[0],1,&ptrans[2*k]);
}

reulermultinom(1,A,&mu_A,1,&adeath);

// Bookkeeping
for (k = 0; k < ESTAGES; k++) {
E[k] -= (etrans[2*k]+etrans[2*k+1]);
E[k+1] += etrans[2*k]; // E[ESTAGES] == L[0]!!
}

E[0] += rpois(b*A);

for (k = 0; k < LSTAGES; k++) {
L[k] -= (ltrans[2*k]+ltrans[2*k+1]);
L[k+1] += ltrans[2*k]; // L[LSTAGES] == P[0]!!
}

for (k = 0; k < PSTAGES; k++) {
P[k] -= (ptrans[2*k]+ptrans[2*k+1]);
P[k+1] += ptrans[2*k]; // P[PSTAGES] == A[0]!!
}
A -= adeath;
"

rmeas_snippet <- 
"
double *E = &E1;
double *L = &L1;
double *P = &P1;

int k;
double E_tot = 0;
double L_tot = 0;
double P_tot = 0;
for (k = 0; k < ESTAGES; k++) E_tot += E[k];
for (k = 0; k < LSTAGES; k++) L_tot += L[k];
for (k = 0; k < PSTAGES; k++) P_tot += P[k];

E_obs = E_tot;
L_obs = L_tot;
P_obs = P_tot;
A_obs = A;"