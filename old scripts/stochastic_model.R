library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
library(magrittr)
library(reshape2)
stopifnot(packageVersion("pomp")>="1.4.9")

# Import and organize the data for the model
read.csv("./data.csv") %>%subset(weeks <= 40, select=c(weeks,rep,L_obs,P_obs,A_obs,DA)) -> dat

dat %>%
  melt(id=c("weeks","rep")) %>%
  acast(variable~rep~weeks) -> datarray

statearray <- datarray
rownames(statearray) <- c("L","P","R","S","A")

# Import and organize the parameters for the model
paramarray <- as.matrix(read.csv('./params.csv'))
row.names(paramarray) <- c("b", "cea", "cel", "cpa", "ua", "ul", "sigma_1", "sigma_2", "sigma_3")
colnames(paramarray) <- 1:24

pomp(
  data = subset(dat, rep==1),
  times="weeks", t0=0,
  initializer=Csnippet("
                       L = 250;
                       P = 5;
                       R = 0;
                       S = 0;
                       A = 100;"),
  rprocess=discrete.time.sim(
    step.fun=Csnippet('
                      R = rbinom(P, exp(-cpa * A));
                      P = rbinom(L, 1 - ul);
                      L = rpois(b*A*exp(-cel * L - cea * A));
                      S = rbinom(A, 1 - ua);
                      A = R + S;'),
    delta.t=2),
  rmeasure=Csnippet("
                    L_obs = L;
                    P_obs = P;
                    A_obs = A;"),
  dmeasure=Csnippet("
                    double eps = 0.000001;
                    if( (abs(L_obs - L) > eps)  ||
                        (abs(P_obs - P) > eps)  ||
                        (abs(A_obs - A) > eps))
                    {
                      lik = 0;
                    } else {
                      lik = 1;
                    }"),
  dprocess=onestep.dens(dens.fun=function(x1,x2,t1,t2,params,...){
    stopifnot(t2==t1+2L)
    with(as.list(params),{
      L1 <- dpois(x2["L"], b*x1["A"]*exp(-cel * x1["L"] - cea * x1["A"]), log=TRUE)
      L2 <- dbinom(x2["P"], x1["L"], 1 - ul, log = TRUE)
      L3 <- dbinom(x2["R"], x1["P"], exp(-cpa * x1["A"]), log = TRUE)
      L4 <- dbinom(x2["S"], x1["A"], 1 - ua, log = TRUE)
      
      L1 + L2 + L3 + L4
    })
    }),
  skeleton=map(
    Csnippet("
             DL = b * A * exp(-cel * L - cea * A);
             DP = L * (1 - ul);
             DA = P * exp(-cpa * A) + A * (1 - ua);"),
    delta.t=2),
  statenames = c("L", "P", "R", "S", "A"),
  paramnames = c("b", "cea", "cel", "cpa", "ua", "ul", "sigma_1", "sigma_2", "sigma_3")) -> model

f2 <- function(par) {
  p <- paramarray
  p[c('b', 'cea','cel','ua','ul','sigma_1','sigma_2','sigma_3'),] <-
    c(par['b'],
      par['cea'],
      par['cel'],
      par['ua'],
      par['ul'],
      par['sigma_1'],
      par['sigma_2'],
      par['sigma_3'])
  sum(dprocess(model,x=statearray,params=p,times=time(model),log=TRUE))
}

optim(fn=f2, control=c(fnscale=-1), par=defaultparams[c('b', 'cea','cel','ua','ul','sigma_1','sigma_2','sigma_3')]) -> fit2
fit2
