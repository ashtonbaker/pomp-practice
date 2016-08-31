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
sum(dprocess(model,x=x,params=defaultparams,times=time(model),log=TRUE))

library(magrittr)
library(reshape2)
raw_data %>% subset(weeks <= 40, select=c(weeks,rep,L,P,A)) %>%
  melt(id=c("weeks","rep")) %>%
  acast(variable~rep~weeks) -> datarray
