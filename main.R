library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
stopifnot(packageVersion("pomp")>="1.4.9")

raw_data <- read.csv("./hunt.csv")
dat = subset(raw_data, rep == 4)
L_0 <- dat[1, "L"]
P_0 <- dat[1, "P"]
A_0 <- dat[1, "A"]

pomp(
  data=dat,
  times="weeks", t0=0,
  skeleton=map(
    Csnippet("
      DL = b * A * exp(-cel * L - cea * A);
      DP = L * (1 - ul);
      DA = P * exp(-cpa * A) + A * (1 - ua);")),
  initializer=Csnippet("
      L = L_0;
      P = P_0;
      A = A_0;"),
  statenames = c("L", "P", "A"),
  paramnames = c("b", "cea", "cel", "cpa", "ua", "ul")) -> deterministic_model
