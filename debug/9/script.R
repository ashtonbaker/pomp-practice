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

# init_snippet <-
#   Csnippet("
#     double *E = &E1;
#     double *L = &L1;
#     double *P = &P1;
#
#     double gamma_E = (ESTAGES / tau_E) * exp((-cel * L_0 - cea * A_0) / ESTAGES);
#     double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
#     double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A_0) / ESTAGES);
#
#     double mu_e = (ESTAGES / tau_E) * (1 - exp((-cel * L_0 - cea * A_0) / ESTAGES));
#     double mu_l = (LSTAGES / tau_L) * mu_L;
#     double mu_p = (PSTAGES / tau_P) * (1 - exp((-cpa * A_0) / ESTAGES));
#
#     double L_rate[LSTAGES] = {0};
#     double P_rate[PSTAGES] = {0};
#
#     int k;
#     double sum;
#     for (k = 0, sum = 0; k < LSTAGES; k++) sum += L_rate[k] = pow(gamma_L/(gamma_L + mu_l), k);
#     for (k = 0; k < LSTAGES; k++) L_rate[k] /= sum;
#     for (k = LSTAGES - 1, sum = 0; k >=0; k--){
#       sum += L_rate[k];
#       L_rate[k] /= sum;
#     }
#
#     for (k = 0, sum = 0; k < PSTAGES; k++) sum += P_rate[k] = pow(gamma_P/(gamma_P + mu_p), k);
#     for (k = 0; k < PSTAGES; k++) P_rate[k] /= sum;
#     for (k = PSTAGES - 1, sum = 0; k >=0; k--){
#       sum += P_rate[k];
#       P_rate[k] /= sum;
#     }
#
#
#     for (k = 0; k < ESTAGES; k++) E[k] = 0;
#
#     int L_count = L_0;
#     for (k = 0; k < LSTAGES - 1; k++){
#       L[k] = rbinom(L_count, L_rate[k]);
#       L_count -= L[k];
#     }
#     L[LSTAGES - 1] = L_count;
#
#     int P_count = P_0;
#     for (k = 0; k < PSTAGES - 1; k++){
#       P[k] = rbinom(P_count, P_rate[k]);
#       P_count -= P[k];
#     }
#     P[PSTAGES - 1] = P_count;
#
#     A = 100;")

init_snippet <-
  Csnippet("
    double *E = &E1;
    double *L = &L1;
    double *P = &P1;
    int k = 0;
    for(k = 0; k < 7; k++) E[k] = 0;
    for(k = 0; k < 5; k++) L[k] = 36;
    for(k = 5; k < 7; k++) L[k] = 35;
    for(k = 0; k < 5; k++) P[k] = 1;
    for(k = 5; k < 7; k++) P[k] = 0;
    A = 100;")

rproc_snippet <-
  Csnippet("
    double *E = &E1;
    double *L = &L1;
    double *P = &P1;

    int k;
    double L_tot = 0;
    for (k = 0; k < LSTAGES; k++) L_tot += L[k];

    double gamma_E = (ESTAGES / tau_E) * exp((-cel * L_tot - cea * A) / ESTAGES);
    double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
    double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A) / ESTAGES);

    double mu_e = (ESTAGES / tau_E) * (1 - exp((-cel * L_tot - cea * A) / ESTAGES));
    double mu_l = (LSTAGES / tau_L) * mu_L;
    double mu_p = (PSTAGES / tau_P) * (1 - exp((-cpa * A) / ESTAGES));

    double etrans[2*ESTAGES], ltrans[2*LSTAGES], ptrans[2*PSTAGES], adeath;

    // Calculate who goes where
    for (k = 0; k < ESTAGES; k++) {
      etrans[2*k]   = rbinom(E[k], gamma_E);                             // Eggs growing to next stage
      etrans[2*k+1] = rbinom( E[k] - etrans[2*k] , mu_e/(1 - gamma_E) ); // Eggs dying
    }

    for (k = 0; k < LSTAGES; k++) {
      ltrans[2*k]   = rbinom( L[k], gamma_L);                          // Larvae growing to next stage
      ltrans[2*k+1] = rbinom( L[k]-ltrans[2*k], mu_l/(1 - gamma_L));   // Larvae dying
    }

    for (k = 0; k < PSTAGES; k++) {
      ptrans[2*k]   = rbinom(P[k], gamma_P);                           // Pupae growing to next stage
      ptrans[2*k+1] = rbinom( P[k]-ptrans[2*k], mu_p/(1 - gamma_P) ); // Pupae dying
    }

    adeath = rbinom(A, mu_A);

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
  const double *L = &L1;
  const double *P = &P1;

  int k;
  double L_tot = 0;
  double P_tot = 0;
  for (k = 0; k < LSTAGES; k++) L_tot += L[k];
  for (k = 0; k < PSTAGES; k++) P_tot += P[k];

  lik = 0;
  lik += dpois(L_obs, L_tot, 1) +
         dpois(P_obs, P_tot, 1) +
         dpois(A_obs, A,     1);

  if(lik == 0){
    printf(\"\\n\\nweeks %f\", t);
    printf(\"\\nL_tot %f\", L_tot);
    printf(\"\\nP_tot %f\", P_tot);
    printf(\"\\nA_tot %f\", A);
    printf(\"\\nL_obs %f\", L_obs);
    printf(\"\\nP_obs %f\", P_obs);
    printf(\"\\nA_obs %f\", A_obs);
  }

  lik = (give_log) ? lik : exp(lik);
   ")

rmeas_snippet <-
  Csnippet("
    const double *E = &E1;
    const double *L = &L1;
    const double *P = &P1;

    int k;
    double L_tot = 0;
    double P_tot = 0;
    for (k = 0; k < LSTAGES; k++) L_tot += L[k];
    for (k = 0; k < PSTAGES; k++) P_tot += P[k];

    L_obs = rpois(L_tot);
    P_obs = rpois(P_tot);
    A_obs = rpois(A);")

from_est <-
  Csnippet("
    Tb = exp(b);
    Tcea = expit(cea);
    Tcel = expit(cel);
    Tcpa = expit(cpa);
    Tmu_A = expit(mu_A);
    Tmu_L = expit(mu_L);
    Ttau_E = exp(tau_E);
    Ttau_L = exp(tau_L);
    Ttau_P = exp(tau_P);")

to_est <-
  Csnippet("
    Tb = log(b);
    Tcea = logit(cea);
    Tcel = logit(cel);
    Tcpa = logit(cpa);
    Tmu_A = logit(mu_A);
    Tmu_L = logit(mu_L);
    Ttau_E = log(tau_E);
    Ttau_L = log(tau_L);
    Ttau_P = log(tau_P);")

pomp(
  data = subset(dat, rep==4),
  times="weeks", t0=0,
  obsnames = c("E_obs", "L_obs", "P_obs", "A_obs"),
  statenames = c(sprintf("E%d",1:stages.E),sprintf("L%d",1:stages.L),sprintf("P%d",1:stages.P),"A"),
  paramnames = c("b", "cea", "cel", "cpa", "mu_A", "mu_L", "tau_E", "tau_L", "tau_P"),
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
             "tau_L"=7.18906255435284,
             "tau_P"=18.0248791283609
             )) -> model

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
                    tau_E=0.02, tau_L=0.02, tau_P=0.02)
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
  tau_E = c(7, 14),
  tau_L = c(7, 14),
  tau_P = c(7, 14),
  tau_A = c(7, 14)
)

print("Starting global search")

stew(file="./output/box_search_global.rda",{
  n_global <- getDoParWorkers()
  t_global <- system.time({
    mf1 <- mifs_local[[1]]
    guesses <- as.data.frame(apply(params_box,1,function(x)runif(300,x[1],x[2])))
    results_global <- foreach(guess=iter(guesses,"row"),
                              .packages='pomp',
                              .combine=rbind,
                              .options.multicore=list(set.seed=TRUE),
                              .export=c("mf1")
    ) %dopar%
    {
      mf <- mif2(mf1,start=c(unlist(guess)))
      mf <- mif2(mf,Nmif=100)
      ll <- replicate(10,logLik(pfilter(mf,Np=100000)))
      ll <- logmeanexp(ll,se=TRUE)
      c(coef(mf),loglik=ll[1],loglik=ll[2])
    }
  })
},seed=1270401374,kind="L'Ecuyer")
results_global <- as.data.frame(results_global)
results <- rbind(results,results_global[names(results)])
write.csv(results,file="./output/model_params.csv",row.names=FALSE)

print("Finished global search")

sink("message.txt", append=FALSE, split=FALSE)
proc.time()
sink()
q(runLast = FALSE)
