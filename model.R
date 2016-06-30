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

statearray <- datarray
rownames(statearray) <- c("L","P","A")

# Import and organize the parameters for the model
paramarray <- as.matrix(read.csv('./params.csv'))
row.names(paramarray) <- c("b", "cea", "cel", "cpa", "ua", "ul", "sigma_1", "sigma_2", "sigma_3")
colnames(paramarray) <- 1:24

params1 <- c(b = 10.45,
                   cea = 0.01310,
                   cel = 0.01731,
                   cpa = 0.004619,
                   ua = 0.007629,
                   ul = 0.2000,
                   sigma_1 = 1.621,
                   sigma_2 = 0.7375,
                   sigma_3 = 0.01212)

params2 <- c(b = 2,
                   cea = 0.5,
                   cel = 0.5,
                   cpa = 0.5,
                   ua = 0.5,
                   ul = 0.5,
                   sigma_1 = 1,
                   sigma_2 = 1,
                   sigma_3 = 1)
pomp(
  data = subset(dat, rep==2),
  times="weeks", t0=0,
  initializer=Csnippet("
                       L = 250;
                       P = 5;
                       A = 100;"),
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
    stopifnot(t2==t1+2L)
    with(as.list(params),{
      mu_l <- sqrt(b * x1["A"] * exp(-cel*x1["L"] - cea*x1["A"]))
      mu_p <- sqrt(x1["L"] * (1 - ul))
      mu_a <- sqrt(x1["P"] * exp(-cpa * x1["A"]) + x1["A"] * (1 - ua))
      
      likl <- dnorm(sqrt(x2["L"]), mean = mu_l, sd = sigma_1,log=TRUE)
      likp <- dnorm(sqrt(x2["P"]), mean = mu_p, sd = sigma_2,log=TRUE)
      lika <- dnorm(sqrt(x2["A"]), mean = mu_a, sd = sigma_3,log=TRUE)
      if(ua < 0){lika = -Inf}
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
  paramnames = c("b", "cea", "cel", "cpa", "ua", "ul", "sigma_1", "sigma_2", "sigma_3")) -> model

f2 <- function(par) {
  p <- paramarray
  p[c('b', 'cea','cel','ul','sigma_1','sigma_2','sigma_3'),] <-
      c(par['b'],
         par['cea'],
         par['cel'],
         par['ul'],
         par['sigma_1'],
         par['sigma_2'],
         par['sigma_3'])
  p['ua',c(4, 11, 24)] <- par['ua']
  p['cpa',c(4, 11, 24)] <- par['cpa']
  sum(dprocess(model,x=statearray,params=p,times=time(model),log=TRUE))
}

optim(fn=f2, control=c(fnscale=-1, maxit=10000), par=params1[c('b', 'cea','cel', 'cpa','ua','ul','sigma_1','sigma_2','sigma_3')]) -> fit2
print(fit2)
