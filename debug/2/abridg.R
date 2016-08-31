library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
library(magrittr)
library(reshape2)
library(foreach)
options(echo = FALSE)
#stopifnot(packageVersion("pomp")>="1.7.4.1")

read.csv("./data/data.csv") %>%subset(weeks <= 40, select=c(weeks,rep,L_obs,P_obs,A_obs)) -> dat
dat$E_obs <- 0
dat %>%
  melt(id=c("weeks","rep")) %>%
  acast(variable~rep~weeks) -> datarray;

stages.E <- 7
stages.L <- 7
stages.P <- 7
stages.A <- 1

glob_snippet <- Csnippet(sprintf("
#include <math.h>
#define ESTAGES %d
#define LSTAGES %d
#define PSTAGES %d
#define ASTAGES %d
#define L_0 250
#define P_0 5
#define A_0 100
", stages.E, stages.L, stages.P, stages.A))

init_snippet <- Csnippet("
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
A = 100;")

rproc_snippet <-
Csnippet("
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

double rate[2], etrans[2*ESTAGES], ltrans[2*LSTAGES], ptrans[2*PSTAGES], adeath;

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
")

dmeas_snippet <-Csnippet(
"
const double *E = &E1;
const double *L = &L1;
const double *P = &P1;

int k;
double E_tot = 0;
double L_tot = 0;
double P_tot = 0;
for (k = 0; k < ESTAGES; k++) E_tot += E[k];
for (k = 0; k < LSTAGES; k++) L_tot += L[k];
for (k = 0; k < PSTAGES; k++) P_tot += P[k];

/*
double eps = 0.000001;
if((abs(L_obs - L_tot) > eps) ||
   (abs(P_obs - P_tot) > eps) ||
   (abs(A_obs - A) > eps)) {
  lik = 0;
} else {
  lik = 1;
}*/

lik =   log(pnorm(L_obs + 0.5, L_tot, meas_sd, 1, 0) - pnorm(L_obs - 0.5, L_tot, meas_sd, 1, 0)) +
        log(pnorm(P_obs + 0.5, P_tot, meas_sd, 1, 0) - pnorm(P_obs - 0.5, P_tot, meas_sd, 1, 0)) +
        log(pnorm(A_obs + 0.5, A,     meas_sd, 1, 0) - pnorm(A_obs - 0.5, A,     meas_sd, 1, 0));

if(isnan(lik))
{
printf(\"\\n\\nL_tot %f\", L_tot);
printf(\"\\nP_tot %f\", P_tot);
printf(\"\\nA_tot %f\", A);
printf(\"\\nsd    %f\", meas_sd);
printf(\"\\nb     %f\", b);
printf(\"\\ncea   %f\", cea);
printf(\"\\ncpa   %f\", cpa);
printf(\"\\nmu_A  %f\", mu_A);
printf(\"\\nmu_L  %f\", mu_L);
printf(\"\\ntau_E %f\", tau_E);
printf(\"\\ntau_L %f\", tau_L);
printf(\"\\ntau_P %f\", tau_P);
}


lik = (give_log) ? lik : exp(lik);

")

rmeas_snippet <-
Csnippet("
const double *E = &E1;
const double *L = &L1;
const double *P = &P1;

int k;
double E_tot = 0;
double L_tot = 0;
double P_tot = 0;
for (k = 0; k < ESTAGES; k++) E_tot += E[k];
for (k = 0; k < LSTAGES; k++) L_tot += L[k];
for (k = 0; k < PSTAGES; k++) P_tot += P[k];

/*
E_obs = E_tot;
L_obs = L_tot;
P_obs = P_tot;
A_obs = A;
*/

E_obs = E_tot;
L_obs = nearbyint(rnorm(L_tot, meas_sd));
P_obs = nearbyint(rnorm(P_tot, meas_sd));
A_obs = nearbyint(rnorm(A, meas_sd));

")

from_est <- Csnippet("
Tb = exp(b);
Tcea = expit(cea);
Tcel = expit(cel);
Tcpa = expit(cpa);
Tmu_A = expit(mu_A);
Tmu_L = expit(mu_L);
Ttau_E = exp(tau_E);
Ttau_L = exp(tau_L);
Ttau_P = exp(tau_P);
Tmeas_sd = exp(meas_sd);
")

to_est <- Csnippet("
Tb = log(b);
Tcea = logit(cea);
Tcel = logit(cel);
Tcpa = logit(cpa);
Tmu_A = logit(mu_A);
Tmu_L = logit(mu_L);
Ttau_E = log(tau_E);
Ttau_L = log(tau_L);
Ttau_P = log(tau_P);
Tmeas_sd = log(meas_sd);
")

pomp(
  data = subset(dat, rep==4),
  times="weeks", t0=0,
  obsnames = c("E_obs", "L_obs", "P_obs", "A_obs"),
  statenames = c(sprintf("E%d",1:stages.E),sprintf("L%d",1:stages.L),sprintf("P%d",1:stages.P),"A"),
  paramnames = c("b", "cea", "cel", "cpa", "mu_A", "mu_L", "tau_E", "tau_L", "tau_P", "meas_sd"),
  globals = glob_snippet,
  initializer = init_snippet,
  rprocess = discrete.time.sim(
    step.fun = rproc_snippet,
    delta.t = 1/7),
  dmeasure = dmeas_snippet,
  rmeasure = rmeas_snippet,
  toEstimationScale = to_est,
  fromEstimationScale = from_est,
  params = c("b"=1.18702207924403,
             "cea"=0.0132088702404268,
             "cel"=0.0172244842038504,
             "cpa"=0.00466955565765198,
             "mu_A"=1.89532307252467e-05,
             "mu_L"=0.0158937470126093,
             "tau_E"=15.7219226675806,
             "tau_L"=5.18906255435284,
             "tau_P"=18.0248791283609,
             "meas_sd" = 10)) -> model

defaultparams <- model@params
pf <- pfilter(model, params = defaultparams, Np=1000)
logLik(pf)

library(foreach);
library(doParallel);
registerDoParallel(cores=40);

print("Starting initial pfilter")

stew(file="./output/pf.rda",{
  t_pf <- system.time(
    pf <- foreach(i=1:10,.packages='pomp',
                  .options.multicore=list(set.seed=TRUE),
                  .export=c("model")
    ) %dopar% {
      pfilter(model,params=defaultparams,Np=10000)
    }
  )
  n_pf <- getDoParWorkers()
},seed=625904618,kind="L'Ecuyer")

print("Finished initial pfilter")

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))
results <- as.data.frame(as.list(c(coef(pf[[1]]),loglik=L_pf[1],loglik=L_pf[2])))
write.csv(results,file="./output/model_params.csv",row.names=FALSE)

print("Starting local box search")

stew(file="./output/box_search_local.rda",{
  t_local_mif <- system.time({
    mifs_local <- foreach(i=1:20,
                          .packages='pomp',
                          .combine=c,
                          .options.multicore=list(set.seed=TRUE),
                          .export=c("model")
    ) %dopar%
    {
      mif2(
        model,
        start=defaultparams,
        Np=2000,
        Nmif=50,
        cooling.type="geometric",
        cooling.fraction.50=0.5,
        transform=TRUE,
        rw.sd=rw.sd(b=0.02, cea=0.02, cel=0.02, cpa=0.02, mu_A=0.02, mu_L=0.02,
                    tau_E=0.02, tau_L=0.02, tau_P=0.02, meas_sd = 0.02)
      )
    }
  })
},seed=482947940,kind="L'Ecuyer")

print("Finished local box search")

print("Starting lik_local")

stew(file="./output/lik_local.rda",{
  t_local_eval <- system.time({
    results_local <- foreach(mf=mifs_local,
                             .packages='pomp',
                             .combine=rbind,
                             .options.multicore=list(set.seed=TRUE)
    ) %dopar%
    {
      evals <- replicate(10, logLik(pfilter(mf,Np=20000)))
      ll <- logmeanexp(evals,se=TRUE)
      c(coef(mf),loglik=ll[1],loglik=ll[2])
    }
  })
},seed=900242057,kind="L'Ecuyer")

print("Finished lik_local")

results_local <- as.data.frame(results_local)
results <- rbind(results,results_local[names(results)])
write.csv(results,file="./output/model_params.csv",row.names=FALSE)

params_box <- rbind(
  b=c(0, 20),
  cea=c(0, 1),
  cel = c(0, 1),
  cpa = c(0, 1),
  mu_A = c(0, 1),
  mu_L = c(0, 1),
  tau_E = c(0, 14),
  tau_L = c(0, 14),
  tau_P = c(0, 14),
  tau_A = c(0, 14),
  meas_sd = c(0, 3)
)

# print("Starting global search")
#
# stew(file="./output/box_search_global.rda",{
#   n_global <- getDoParWorkers()
#   t_global <- system.time({
#     mf1 <- mifs_local[[1]]
#     guesses <- as.data.frame(apply(params_box,1,function(x)runif(300,x[1],x[2])))
#     results_global <- foreach(guess=iter(guesses,"row"),
#                               .packages='pomp',
#                               .combine=rbind,
#                               .options.multicore=list(set.seed=TRUE),
#                               .export=c("mf1")
#     ) %dopar%
#     {
#       mf <- mif2(mf1,start=c(unlist(guess)))
#       mf <- mif2(mf,Nmif=100)
#       ll <- replicate(10,logLik(pfilter(mf,Np=100000)))
#       ll <- logmeanexp(ll,se=TRUE)
#       c(coef(mf),loglik=ll[1],loglik=ll[2])
#     }
#   })
# },seed=1270401374,kind="L'Ecuyer")
# results_global <- as.data.frame(results_global)
# results <- rbind(results,results_global[names(results)])
# write.csv(results,file="./output/model_params.csv",row.names=FALSE)
#
# print("Finished global search")

sink("message.txt", append=FALSE, split=FALSE)
proc.time()
sink()
q(runLast = FALSE)
