library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
library(magrittr)
library(reshape2)
stopifnot(packageVersion("pomp")>="1.4.9")

# Import and organize the data for the model
read.csv("./data.csv") %>%subset(weeks <= 40, select=c(weeks,rep,L_obs,P_obs,A_obs)) -> dat

dat %>%
  melt(id=c("weeks","rep")) %>%
  acast(variable~rep~weeks) -> datarray

pomp(
  data = subset(dat, rep==2),
  times="weeks", t0=0,
  globals = "int lstages = 14, pstages = 14, astages = 1, L_0 = 250, P_0 = 5, A_0 = 100;",
  initializer=Csnippet("
                      double mu_P = cpa * A;
                      double *L = &L1;
                      double *P = &P1;
                      double L_rate[100], P_rate[100] = {0};
                      
                      int k;
                      for (k = 0; k < lstages; k++) L_rate[k] = pow(gamma_L / (gamma_L + mu_L), k);
                      for (k = 0; k < pstages; k++) P_rate[k] = pow(gamma_P / (gamma_P + mu_P), k);

                      reulermultinom(lstages, L_0, &L_rate[0], 1, &L[0]);
                      reulermultinom(pstages, P_0, &P_rate[0], 1, &P[0]);
                      A = 100;"),
  statenames = c(sprintf("L%d",1:14),sprintf("P%d",1:14),"A"),
  rprocess=discrete.time.sim(
    step.fun=Csnippet('
                      double *L = &L1;
                      double *P = &P1;
                      double rate[2], ltrans[28], ptrans[28], adeath;
                      int k;
                      double eggs = rpois(b*A*dt/2);
                      for (k = 0; k < lstages; k++) {
                        rate[0] = gamma_L*lstages;
                        rate[1] = mu_L;
                        reulermultinom(2,L[k],&rate[0],dt/2,&ltrans[2*k]);
                      }
                      for (k = 0; k < pstages; k++) {
                        rate[0] = gamma_P*pstages;
                        rate[1] = cpa*A;
                        reulermultinom(2,P[k],&rate[0],dt/2,&ptrans[2*k]);
                      }
                      reulermultinom(1,A,&mu_A,dt/2,&adeath);
                      L[0] += eggs;
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
                    const double *L = &L1;
                    const double *P = &P1;
                    int k;
                    int L_tot = 0, P_tot = 0;
                    for (k = 0; k < lstages; k++) L_tot += L[k];
                    for (k = 0; k < pstages; k++) P_tot += P[k];

                    L_obs = L_tot;
                    P_obs = P_tot;
                    A_obs = A;"),
  dmeasure=Csnippet("
                    const double *L = &L1;
                    const double *P = &P1;
                    int k, L_tot, P_tot;
                    for (k = 0; L_tot = 0; k < lstages, k++) L_tot += L[k];
                    for (k = 0; P_tot = 0; k < pstages, k++) P_tot += P[k];

                    double eps = 0.000001;
                    if((abs(L_obs - L_tot) > eps) ||
                    (abs(P_obs - P_tot) > eps) ||
                    (abs(A_obs - A) > eps)) {
                    lik = 0;
                    } else {
                    lik = 1;
                    }"),
  paramnames = c("b", "cea", "cel", "cpa", "mu_A", "mu_L", "gamma_L", "gamma_P", "sigma_1", "sigma_2", "sigma_3"),
  params = c("b"=10.45, "cea"=0.01310, "cel"=0.01731, "cpa"=0.004619, "mu_A"=0.007629,"mu_L"=0.2,"gamma_L"=0.1,"gamma_P"=0.1,"sigma_1"=1.621, "sigma_2"=0.7375, "sigma_3"=0.01212)) -> model

