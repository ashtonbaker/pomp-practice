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

source("./snippets.R")

pomp(
  data = subset(dat, rep==4),
  times="weeks", t0=0,
  obsnames = c("E_obs", "L_obs", "P_obs", "A_obs"),
  statenames = c(sprintf("E%d",1:5),sprintf("L%d",1:7),sprintf("P%d",1:7),"A"),
  paramnames = c("b", "cea", "cel", "cpa", "mu_A", "mu_L", "tau_E", "tau_L", "tau_P", "tau_A"),
  globals = glob_snippet,
  initializer=Csnippet(init_snippet),
  rprocess=discrete.time.sim(
    step.fun=Csnippet(rproc_snippet),
    delta.t=1/7),
  rmeasure=Csnippet(rmeas_snippet),
  params = c("b"=0.7464286,
             "cea"=0.01310,
             "cel"=0.01731,
             "cpa"=0.004619,
             "mu_A"=0.007629,
             "mu_L"=0.015812,
             "tau_E"=5,
             "tau_L"=14,
             "tau_P"=14,
             "tau_A"=14)) -> model

defaultparams <- model@params
  
ssr <- function(par) {
  total <- 0
  if(min(defaultparams) < 0){
    total <- -Inf
  }
  
  if(max(c(cea, cel, cpa) > 1){
    total <- -Inf
  }

  sim <- simulate(model, nsim = 200, params = par)

  for(i in 1:200){
    for(j in c(4)){ #, 11, 24)){
      total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
    }
  }
  
  # a = par
  # a['cpa'] <- 0.00
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(5, 12, 15)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  # 
  # a['cpa'] <- 0.05
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(1, 7, 20)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  # 
  # a['cpa'] <- 0.10
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(6, 10, 16)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  # 
  # a['cpa'] <- 0.25
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(8, 17, 21)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  # 
  # a['cpa'] <- 0.35
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(2, 13, 22)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  # 
  # a['cpa'] <- 0.50
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(9, 14, 19)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  # 
  # a = par
  # a['cpa'] <- 1.00
  # a['mu_A'] <- 0.96
  # sim <- simulate(model, nsim = 200, params = a)
  # for(i in 1:200){
  #   for(j in c(3, 18, 23)){
  #     total = total + sum(rowSums((datarray[c("L_obs", "P_obs", "A_obs"),j,] - obs(sim[[i]])[c('L_obs','P_obs','A_obs'),])^2)^0.5)
  #   }
  # }
  total
}

fit1 <- optim(fn=ssr, par=defaultparams, control=c(maxit=100000))
model <- pomp(model, params=fit1$par)

sims <- simulate(model, nsim=40,as.data.frame=TRUE,include.data=TRUE)

ggplot(sims,mapping=aes(x=time,y=L_obs,group=sim,color=sim=="data"))+geom_line()+guides(color=FALSE)
ggplot(sims,mapping=aes(x=time,y=P_obs,group=sim,color=sim=="data"))+geom_line()+guides(color=FALSE)
ggplot(sims,mapping=aes(x=time,y=A_obs,group=sim,color=sim=="data"))+geom_line()+guides(color=FALSE)
