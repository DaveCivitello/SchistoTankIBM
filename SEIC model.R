# SEIC model for Fig 1 in mesocosm paper

setwd("C:/RData")
library(deSolve)

SEIC <- function(t,y,params){
 S = y[1]; E = y[2]; I = y[3]; C = y[4];
  with(
    as.list(params),
    {
      dS <-  b*(1 - (S+E+I)/K)*(S+E) - mu_s*S - BM*S
      dE <- BM*S - mu_s*E - sigma*E
      dI <- sigma*E - (mu_s + mu_i)*I
      dC <-lambda*I - mu_c*C
      res <- c(dS, dE, dI, dC)
      list(res)
    }
  )
}

initials = c(S = 0.1, E = 0, I = 0, C = 0)
parameters = c(b = 0.1, K = 5, mu_s = 0.01, BM = 0.01, sigma = 1/28, mu_i = 0.04, lambda = 50, mu_c = 1)

sim1 = data.frame(ode(y = initials, times = 0:100, parms = parameters, func = SEIC))
sim1[,"N"] <- sim1[,"S"] + sim1[,"E"] + sim1[,"I"]
plot(N ~ time, data=sim1, type="l", ylim=c(0, 20))
lines(I ~ time, data=sim1, col="blue")
lines(C ~ time, data=sim1, col="red")

