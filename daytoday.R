library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
library(magrittr)
library(reshape2)
stopifnot(packageVersion("pomp")>="1.4.9")

# Import and organize the data for the model
read.csv("./data.csv") %>%subset(weeks <= 40, select=c(weeks,rep,L_obs,P_obs,A_obs)) -> dat
dat$E_obs <- 0


dat %>%
  melt(id=c("weeks","rep")) %>%
  acast(variable~rep~weeks) -> datarray

pomp(
  data = subset(dat, rep==2),
  times="weeks", t0=0,
  obsnames = c("E_obs", "L_obs", "P_obs", "A_obs"),
  statenames = c(sprintf("E%d",1:14),sprintf("L%d",1:14),sprintf("P%d",1:14),"A"),
  paramnames = c("b", "cea", "cel", "cpa", "mu_A", "mu_L", "gamma_L", "sigma_1", "sigma_2", "sigma_3"),
  globals = "int estages = 14, lstages = 14, pstages = 14, astages = 1, L_0 = 250, P_0 = 5, A_0 = 100;",
  initializer=Csnippet("
                      #define LSTAGE 14
                      double *E = &E1;
                      double *L = &L1;
                      double *P = &P1;
                      double L_rate[LSTAGE], P_rate[100] = {0};

                      int k;
                      double L_tot = 0;
                      for (k = 0; k < lstages; k++) L_tot += L[k];
                      
                      double gamma_E = exp((-0.01731 * L_tot - 0.01310 * A) * 0.07142857);
                      double mu_E = 1 - gamma_E; //pow((1 - exp(-0.01731 * L_tot - 0.01310 * A)),0.07142857);
                      double gamma_P = exp(-0.004619 * A / 14.0);
                      double mu_P = 1 - gamma_P;//pow((1 - exp(-0.004619 * A)),1.0/14.0);
                      
                      for (k = 0; k < lstages; k++) L_rate[k] = pow(gamma_L, k);
                      for (k = 0; k < pstages; k++) P_rate[k] = pow(gamma_P, k);

                      for (k = 0; k < estages; k++) E[k] = 0;
                      reulermultinom(lstages, L_0, &L_rate[0], 1, &L[0]);
                      reulermultinom(pstages, P_0, &P_rate[0], 1, &P[0]);
                      A = 100;"),
  rprocess=discrete.time.sim(
    step.fun=Csnippet('
                      double *E = &E1;
                      double *L = &L1;
                      double *P = &P1;
                      double rate[2], etrans[28], ltrans[28], ptrans[28], adeath;
                      int k;
                      double L_tot = 0;
                      for (k = 0; k < lstages; k++) L_tot += L[k];
                      double gamma_E = exp((-0.01731 * L_tot - 0.01310 * A) * 0.07142857);
                      double mu_E = 1 - gamma_E;//pow((1 - exp(-0.01731 * L_tot - 0.01310 * A)),0.07142857);
                      double gamma_P = exp(-0.004619 * A / 14.0);
                      double mu_P = 1 -gamma_P;//pow((1 - exp(-0.004619 * A)),1.0/14.0);
                      
                      // Calculate who goes where
                      for (k = 0; k < estages; k++) {
                        rate[0] = gamma_E;
                        rate[1] = mu_E;
                        reulermultinom(2,E[k],&rate[0],1,&etrans[2*k]);
                      }

                      for (k = 0; k < lstages; k++) {
                        rate[0] = gamma_L;
                        rate[1] = mu_L;
                        reulermultinom(2,L[k],&rate[0],1,&ltrans[2*k]);
                      }

                      for (k = 0; k < pstages; k++) {
                        rate[0] = gamma_P;
                        rate[1] = mu_P;
                        reulermultinom(2,P[k],&rate[0],1,&ptrans[2*k]);
                      }

                      reulermultinom(1,A,&mu_A,1,&adeath);

                      // Bookkeeping
                      for (k = 0; k < estages; k++) {
                        E[k] -= (etrans[2*k]+etrans[2*k+1]);
                        E[k+1] += etrans[2*k]; // E[estages] == L[0]!!
                      }
                      
                      E[0] += rpois(b*A);

                      for (k = 0; k < lstages; k++) {
                        L[k] -= (ltrans[2*k]+ltrans[2*k+1]);
                        L[k+1] += ltrans[2*k]; // L[lstages] == P[0]!!
                      }

                      for (k = 0; k < pstages; k++) {
                        P[k] -= (ptrans[2*k]+ptrans[2*k+1]);
                        P[k+1] += ptrans[2*k]; // P[pstages] == A[0]!!
                      }
                      A -= adeath;
                     '),
    delta.t=1/7),
  rmeasure=Csnippet("
                    const double *E = &E1;
                    const double *L = &L1;
                    const double *P = &P1;
                    int k;
                    double E_tot = 0, L_tot = 0, P_tot = 0;
                    for (k = 0; k < estages; k++){ E_tot += E[k]; }
                    for (k = 0; k < lstages; k++){ L_tot += L[k]; }
                    for (k = 0; k < pstages; k++){ P_tot += P[k]; }
                    
                    E_obs = E_tot;
                    L_obs = L_tot;
                    P_obs = P_tot;
                    A_obs = A;"),
  params = c("b"=0.7464286,
             "cea"=0.01310,
             "cel"=0.01731,
             "cpa"=0.004619,
             "mu_A"=0.007629,
             "mu_L"=0.015812,
             "gamma_L"=0.984188,
             "sigma_1"=1.621,
             "sigma_2"=0.7375,
             "sigma_3"=0.01212)) -> model

defaultparams <- c("b"=0.7464286, "cea"=0.01310, "cel"=0.01731, "cpa"=0.004619, "mu_A"=0.007629, "mu_L"=0.015812, "gamma_L"=0.984188, "sigma_1"=1.621, "sigma_2"=0.7375, "sigma_3"=0.01212)

ssr <- function(par) {
  total <- 0
  if(sum(defaultparams < 0) != 0){
    total <- -Inf
  }
  
  sim <- simulate(model, nsim = 200, params = par)

  for(i in 1:200){
    for(j in c(4, 11, 24)){
      total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
    }
  }
  total
}

f1 <- function(b_hat) {
  p <- defaultparams
  p['b'] <- b_hat
  ssr(p)
}

fit1 <- optim(fn=ssr, par=defaultparams)