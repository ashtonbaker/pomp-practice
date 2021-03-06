library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
stopifnot(packageVersion("pomp")>="1.4.9")

raw_data <- read.csv("./hunt.csv")
dat = subset(raw_data, rep == 4, select=c(weeks, L, P, A))
names(dat) <- c('weeks', 'L_obs', 'P_obs', "A_obs")

defaultparams <- c(b = 10.45,
                   cea = 0.01310,
                   cel = 0.01731,
                   cpa = 0.004619,
                   ua = 0.007629,
                   ul = 0.2000,
                   L_0 = dat[1, "L_obs"],
                   P_0 = dat[1, "P_obs"],
                   A_0 = dat[1, "A_obs"],
                   sigma_1 = 1.621,
                   sigma_2 = 0.7375,
                   sigma_3 = 0.01212)

mle_params <- c(b = 11.1,
                cea = 0.013,
                cel = 0.014,
                cpa = 0.004,
                ua = 0.003,
                ul = 0.312,
                L_0 = 250,
                P_0 = 5,
                A_0 = 100,
                sigma_1 = 1.325,
                sigma_2 = 0.770,
                sigma_3 = 0.107)

pomp(
  data=dat,
  times="weeks", t0=-1,
  initializer=Csnippet("
                       L = L_0;
                       P = P_0;
                       A = A_0;"),
  rprocess=discrete.time.sim(
    step.fun=Csnippet('
                      double e1 = rnorm(0,sigma_1);
                      double e2 = rnorm(0,sigma_2);
                      double e3 = rnorm(0,sigma_3);
                      L = (sqrt(b * A * exp(-cel * L - cea * A)) + e1)*(sqrt(b * A * exp(-cel * L - cea * A)) + e1);
                      P = (sqrt(L * (1 - ul)) + e2)*(sqrt(L * (1 - ul)) + e2);
                      A = (sqrt(P * exp(-cpa * A) + A * (1 - ua)) + e3)*(sqrt(P * exp(-cpa * A) + A * (1 - ua)) + e3);'),
    delta.t=2),
  rmeasure=Csnippet("
                    L_obs = L;
                    P_obs = P;
                    A_obs = A;"),
  dmeasure=Csnippet("
                    double eps = 0.000001;
                    if((abs(L_obs - L) > eps) ||
                    (abs(P_obs - P) > eps) ||
                    (abs(A_obs - A) > eps)) {
                    lik = 0;
                    } else {
                    lik = 1;
                    }"),
  dprocess=onestep.dens(dens.fun=function(x1,x2,t1,t2,params,...){
    #stopifnot(t2==t1+2L)
    with(as.list(params),{
      mu_l <- sqrt(b * x1["A"] * exp(-cel*x1["L"] - cea*x1["A"]))
      mu_p <- sqrt(x1["L"] * (1 - ul))
      mu_a <- sqrt(x1["P"] * exp(-cpa * x1["A"]) + x1["A"] * (1 - ua))
      
      likl <- dnorm(sqrt(x2["L"]), mean = mu_l, sd = sigma_1,log=TRUE)
      likp <- dnorm(sqrt(x2["P"]), mean = mu_p, sd = sigma_2,log=TRUE)
      lika <- dnorm(sqrt(x2["A"]), mean = mu_a, sd = sigma_3,log=TRUE)
      likl + likp + lika
    })
    }),
  skeleton=map(
    Csnippet("
             DL = b * A * exp(-cel * L - cea * A);
             DP = L * (1 - ul);
             DA = P * exp(-cpa * A) + A * (1 - ua);"),
    delta.t=2),
  statenames = c("L", "P", "A"),
  paramnames = c("b", "cea", "cel", "cpa", "ua", "ul", "L_0", "P_0", "A_0", "sigma_1", "sigma_2", "sigma_3")) -> model


obs(model) -> x
rownames(x) <- c("L","P","A")


dprocess(model,x=x,params=defaultparams,times=time(model),log=TRUE)
xx <- simulate(model,params=defaultparams)
dprocess(xx,x=states(xx),params=coef(xx),times=time(xx),log=TRUE)
## log likelihood of the full set of state-transitions:

#f1 <- function(sigma_3) {
#  p <- defaultparams
#  p['sigma_3'] <- sigma_3
#  sum(dprocess(model,x=x,params=p,times=time(model),log=TRUE))
#}

#sigma_3 <- seq(from=0, to=2, by=0.001)
#LIK <- sapply(sigma_3, f1)
#sigma_3.hat <- sigma_3[which.max(LIK)]
#plot(sigma_3, LIK, type='l')
#abline(v=sigma_3.hat, lty=2)

f2 <- function(par) {
  p <- c(b = par[1],
         cea = par[2],
         cel = par[3],
         cpa = par[4],
         ua = par[5],
         ul = par[6],
         L_0 = 250,
         P_0 = 5,
         A_0 = 100,
         sigma_1 = par[7],
         sigma_2 = par[8],
         sigma_3 = par[9])
  -sum(dprocess(model,x=x,params=p,times=time(model),log=TRUE))
}

optim(fn=f2, par=c(b.hat, cea.hat, cel.hat, cpa.hat, ua.hat,ul.hat,sigma_1.hat,sigma_2.hat,sigma_3.hat)mle) -> fit2
fit2