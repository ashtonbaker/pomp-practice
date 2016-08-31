library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
stopifnot(packageVersion("pomp")>="1.4.9")

raw_data <- read.csv("./hunt.csv")
dat = subset(raw_data, rep == 4, select=c(weeks, L, P, A))
names(dat) <- c('weeks', 'Larvae', 'Pupae', "Adults")

defaultparams <- c(b = 10.45,
                   cea = 0.01310,
                   cel = 0.01731,
                   cpa = 0.004619,
                   ua = 0.007629,
                   ul = 0.2000,
                   L_0 = dat[1, "Larvae"],
                   P_0 = dat[1, "Pupae"],
                   A_0 = dat[1, "Adults"])

pomp(
  data=dat,
  times="weeks", t0=0,
  skeleton=map(
    Csnippet("
      DL = b * A * exp(-cel * L - cea * A);
      DP = L * (1 - ul);
      DA = P * exp(-cpa * A) + A * (1 - ua);"),delta.t=2),
  initializer=Csnippet("
      L = L_0;
      P = P_0;
      A = A_0;"),
  statenames = c("L", "P", "A"),
  paramnames = c("b", "cea", "cel", "cpa", "ua", "ul", "L_0", "P_0", "A_0")) -> deterministic_model

x <- trajectory(deterministic_model, params=defaultparams, as.data.frame=TRUE)
#plot(A~time, data=x,type='o')

ggplot(data=join(as.data.frame(deterministic_model),x,by='time'),
       mapping=aes(x=time))+
  geom_line(aes(y=Larvae),color='black')+
  geom_line(aes(y=L),color='red')

ggplot(data=join(as.data.frame(deterministic_model),x,by='time'),
       mapping=aes(x=time))+
  geom_line(aes(y=Pupae),color='black')+
  geom_line(aes(y=P),color='red')

print(ggplot(data=join(as.data.frame(deterministic_model),x,by='time'),
       mapping=aes(x=time))+
  geom_line(aes(y=Adults),color='black')+
  geom_line(aes(y=A),color='red')
)
