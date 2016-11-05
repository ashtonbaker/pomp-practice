library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
library(magrittr)
library(reshape2)
library(foreach)
#options(echo = FALSE)
library(doMPI)
library(doRNG)

cl <- startMPIcluster(maxcores = 40,workdir = '~/pomp-practice/debug/naive-looping/')
registerDoMPI(cl)

setwd('~/pomp-practice/debug/naive-looping/')

optsN <- list(123, normal.kind="Ahrens")

stopifnot(packageVersion("pomp")>="1.8.8.1")

read.csv("./data/data.csv") %>%
  subset(weeks <= 40, select=c(weeks,rep,L_obs,P_obs,A_obs)) -> dat
dat %>%
  melt(id=c("weeks","rep")) %>%
  acast(variable~rep~weeks) -> datarray

stages.E <- 1
stages.L <- 6
stages.P <- 7
stages.A <- 1

cpa = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
cpa[c(4, 11, 24)] = 0
cpa[c(5, 12, 15)] = 0.00
cpa[c(1, 7, 20)] = 0.05
cpa[c(6, 10, 16)] = 0.10
cpa[c(8, 17, 21)] = 0.25
cpa[c(2, 13, 22)] = 0.35
cpa[c(9, 14, 19)] = 0.50
cpa[c(3, 18, 23)] = 1.00

mu_A = c(0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96)
mu_A[c(4, 11, 24)] = 0.0

for (i in 1:24) {
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
      A = A_0;
      A_prev = A_0;
      P_prev = P_0;")

  rproc_snippet <-
    Csnippet(sprintf("
      double *E = &E1;
      double *L = &L1;
      double *P = &P1;

      int k;
      double L_tot = 0;
      for (k = 0; k < LSTAGES; k++) L_tot += L[k];

      double gamma_E = (ESTAGES / tau_E) *
                       exp((-cel * L_tot - cea * A) / ESTAGES);
      double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
      double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A) / PSTAGES);

      double mu_e = (ESTAGES / tau_E) - gamma_E;
      double mu_l = (LSTAGES / tau_L) - gamma_L;
      double mu_p = (PSTAGES / tau_P) - gamma_P;

      double etrans[2*ESTAGES], ltrans[2*LSTAGES], ptrans[2*PSTAGES], adeath;

      // Calculate who goes where
      for (k = 0; k < ESTAGES; k++) {
        // Eggs growing to next stage
        etrans[2*k]   = rbinom(E[k], gamma_E);

        // Eggs dying
        etrans[2*k+1] = rbinom(E[k]-etrans[2*k], mu_e/(1 - gamma_E) );
      }

      for (k = 0; k < LSTAGES; k++) {
        // Larvae growing to next stage
        ltrans[2*k]   = rbinom(L[k], gamma_L);

        // Larvae dying
        ltrans[2*k+1] = rbinom(L[k]-ltrans[2*k], mu_l/(1 - gamma_L));
      }

      for (k = 0; k < PSTAGES; k++) {
        // Pupae growing to next stage
        ptrans[2*k]   = rbinom(P[k], gamma_P);

        // Pupae dying
        ptrans[2*k+1] = rbinom(P[k]-ptrans[2*k], mu_p/(1 - gamma_P) );
      }

      adeath = rbinom(A, mu_A);

      // Bookkeeping
      E[0] += rpois(b*A); // oviposition

      for (k = 0; k < ESTAGES; k++) {
        // Subtract eggs that die or progress
        E[k] -= (etrans[2*k]+etrans[2*k+1]);

        // Add eggs that arrive from previous E stage.
        E[k+1] += etrans[2*k]; // E[ESTAGES] == L[0]!!
      }

      for (k = 0; k < LSTAGES; k++) {
        // Subtract larvae that die or progress
        L[k] -= (ltrans[2*k]+ltrans[2*k+1]);

        // Add larvae that arrive from previous E stage.
        L[k+1] += ltrans[2*k]; // L[LSTAGES] == P[0]!!
      }

      for (k = 0; k < PSTAGES; k++) {
        // Subtract pupae that die or progress
        P[k] -= (ptrans[2*k]+ptrans[2*k+1]);

        // Add pupae that arrive from previous E stage.
        P[k+1] += ptrans[2*k]; // P[PSTAGES] == A[0]!!
      }

      A -= adeath;

      if ((nearbyint(t + 0.01) %% 14 == 0) && (t != 0) && %f > 0.5) {
        double P_tot = 0;
        for (k = 0; k < PSTAGES; k++) P_tot += P[k];

        double A_pred = nearbyint((1 - 0.96) * A_prev) + nearbyint(P_prev * exp(-%f * A));
        if (A_pred < A) {
          A = A_pred;
        }
        P_prev = P_tot;
        A_prev = A;
      }
      ", mu_A[i], cpa[i]))

  dmeas_snippet <- Csnippet(
    "
    const double *L = &L1;
    const double *P = &P1;
    double fudge = 1e-9;

    int k;
    double L_tot = 0;
    double P_tot = 0;
    for (k = 0; k < LSTAGES; k++) L_tot += L[k];
    for (k = 0; k < PSTAGES; k++) P_tot += P[k];

    lik = dnbinom_mu(L_obs, 1/od, L_tot+fudge, 1) +
          dnbinom_mu(P_obs, 1/od, P_tot+fudge, 1) +
          dnbinom_mu(A_obs, 1/od, A+fudge,     1);

  //  if(lik < -138){
  //    Rprintf(\"\\n\\nweeks %f\", t);
  //    Rprintf(\"\\nL_tot %f\", L_tot);
  //    Rprintf(\"\\nP_tot %f\", P_tot);
  //    Rprintf(\"\\nA_tot %f\", A);
  //    Rprintf(\"\\nL_obs %f\", L_obs);
  //    Rprintf(\"\\nP_obs %f\", P_obs);
  //    Rprintf(\"\\nA_obs %f\", A_obs);
  //    Rprintf(\"\\nloglik %f\",lik);
  //  }

    lik = (give_log) ? lik : exp(lik);
     ")

  rmeas_snippet <-
    Csnippet("
      const double *L = &L1;
      const double *P = &P1;
      double fudge = 1e-9;

      int k;
      double L_tot = 0;
      double P_tot = 0;
      for (k = 0; k < LSTAGES; k++) L_tot += L[k];
      for (k = 0; k < PSTAGES; k++) P_tot += P[k];

      L_obs = rnbinom_mu(1/od,L_tot+fudge);
      P_obs = rnbinom_mu(1/od,P_tot+fudge);
      A_obs = rnbinom_mu(1/od,A+fudge);")

  from_est <-
    Csnippet("
      Tb = exp(b);
      Tcea = expit(cea);
      Tcel = expit(cel);
      Tcpa = expit(cpa);
      Tmu_A = expit(mu_A);
      Tmu_L = expit(mu_L);
      Ttau_E = ESTAGES+exp(tau_E);
      Ttau_L = LSTAGES+exp(tau_L);
      Ttau_P = PSTAGES+exp(tau_P);
      Tod = exp(od);")

  to_est <-
    Csnippet("
      Tb = log(b);
      Tcea = logit(cea);
      Tcel = logit(cel);
      Tcpa = logit(cpa);
      Tmu_A = logit(mu_A);
      Tmu_L = logit(mu_L);
      Ttau_E = log(tau_E-ESTAGES);
      Ttau_L = log(tau_L-LSTAGES);
      Ttau_P = log(tau_P-PSTAGES);
      Tod = log(od);")

  pomp(
    data = subset(dat, rep==i, select=-rep),
    times="weeks", t0=0,
    statenames = c(sprintf("E%d",1:stages.E),
                   sprintf("L%d",1:stages.L),
                   sprintf("P%d",1:stages.P),
                   "A",
                   "A_prev",
                   "P_prev"),
    paramnames = c("b", "cea", "cel", "cpa", "mu_A", "mu_L",
                   "tau_E", "tau_L", "tau_P","od"),
    globals = glob_snippet,
    initializer = init_snippet,
    rprocess = discrete.time.sim(
      step.fun = rproc_snippet,
      delta.t = 1/7),
    dmeasure = dmeas_snippet,
    rmeasure = rmeas_snippet,
    toEstimationScale = to_est,
    fromEstimationScale = from_est,
    params = c(b=1.18702207924403,
               cea=0.0132088702404268,
               cel=0.0172244842038504,
               cpa=0.00466955565765198,
               mu_A=1.89532307252467e-05,
               mu_L=0.0158937470126093,
               tau_E=15.7219226675806,
               tau_L=7.18906255435284,
               tau_P=18.0248791283609,
               od = 1
               )) -> model

  #model %>% simulate(as.data.frame=T,nsim=5) %>%
  #  melt(id=c("time","sim")) %>%
  #  subset(variable %in% c("L_obs","P_obs","A_obs")) %>%
  #  ggplot(aes(x=time,y=value,color=variable,group=sim))+
  #    geom_line()+
  #    facet_wrap(~variable,ncol=1,scales="free_y")

  #pf <- pfilter(model, Np=100000)
  #logLik(pf)


  print("Starting initial pfilter")

  stew(file=sprintf("./output/pf%d.rda", i),{
    t_pf <- system.time(
      pf <- foreach(i=1:10,
                    .packages='pomp',
                    .options.RNG = optsN,
                    .export=c("model")
      ) %dorng% {
        pfilter(model,Np=1000)
      }
    )
    n_pf <- getDoParWorkers()
  },seed=625904618,kind="L'Ecuyer")

  print("Finished initial pfilter")

  (L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))
  results <- as.data.frame(as.list(c(coef(pf[[1]]),loglik=L_pf[1],loglik=L_pf[2])))
  write.csv(results,file="./output/model_params.csv",row.names=FALSE)

  print("Starting local box search")

  stew(file=sprintf("./output/box_search_local%i.rda", i),{
    t_local_mif <- system.time({
      mifs_local <- foreach(i=1:20,
                            .packages='pomp',
                            .options.RNG = optsN,
                            .combine=c,
                            .export=c("model")
      ) %dorng%
      {
        mif2(
          model,
          Np=1000,
          Nmif=20,
          cooling.type="geometric",
          cooling.fraction.50=0.5,
          transform=TRUE,
          rw.sd=rw.sd(b=0.02, cea=0.02, cel=0.02, cpa=0.02,
                      mu_A=0.02, mu_L=0.02, od=0.02,
                      tau_E=0.02, tau_L=0.02, tau_P=0.02)
        )
      }
    })
  },seed=482947940,kind="L'Ecuyer")

  print("Finished local box search")

  print("Starting lik_local")

  stew(file=sprintf("./output/lik_local%d.rda", i),{
    t_local_eval <- system.time({
      results_local <- foreach(mf=mifs_local,
                               .options.RNG = optsN,
                               .packages='pomp',
                               .combine=rbind
      ) %dorng%
      {
        evals <- replicate(10, logLik(pfilter(mf,Np=1000)))
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
    b=c(0, 2),
    cea=c(0, 1),
    cel = c(0, 1),
    cpa = c(0, 1),
    mu_A = c(0, 1),
    mu_L = c(0, 1),
    tau_E = c(7, 14),
    tau_L = c(7, 14),
    tau_P = c(7, 14),
    od = c(1,1)
  )

  print("Starting global search")

  stew(file=sprintf("./output/box_search_global%d.rda", i),{
    n_global <- getDoParWorkers()
    t_global <- system.time({
      mf1 <- mifs_local[[1]]
      guesses <- as.data.frame(apply(params_box,1,function(x)runif(30,x[1],x[2])))
      results_global <- foreach(guess=iter(guesses,"row"),
                                .options.RNG = optsN,
                                .packages='pomp',
                                .combine=rbind,
                                .export=c("mf1")
      ) %dorng%
      {
        mf <- mif2(mf1,start=c(unlist(guess)),tol=1e-60)
        mf <- mif2(mf,Nmif=100)
        ll <- replicate(10,logLik(pfilter(mf,Np=10000)))
        ll <- logmeanexp(ll,se=TRUE)
        c(coef(mf),loglik=ll[1],loglik=ll[2])
      }
    })
  },seed=1270401374,kind="L'Ecuyer")
  results_global <- as.data.frame(results_global)
  results <- rbind(results,results_global[names(results)])
  #write.csv(results,file="./output/model_params.csv",row.names=FALSE)

  print("Finished global search")

  p_optim <- results_global[which.max(results_global$loglik),]
  print(p_optim)
  write.table(p_optim, file = "./output/optim_params.csv", append = TRUE, col.names=(i==1), row.names = FALSE, sep=", ")
}

closeCluster(cl)
mpi.quit()
