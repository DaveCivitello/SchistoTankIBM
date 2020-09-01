
# Loading packages
library(Matrix)
library(deSolve)
library("mvtnorm")
library(LaplacesDemon)
library(coda)
library(adaptMCMC)

# Rtools is needed for compiling the .c file of the DEB model
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

# Working directory (put .c file and parameters file here)
setwd("C:/RData")

# compile my model from C definition
dyn.unload("IndividualModel_IBM3.dll") # unload dll
system("R CMD SHLIB IndividualModel_IBM3.c")
dyn.load("IndividualModel_IBM3.dll") # Load dll

# DEB Parameter estimates current as of Civitello et al. 2020 Proc B (Highest posterior density parameters)
pars = readRDS("Starvation_parameters.RDa")

# Adds parameters for the environment, resource production model, and transmission model
pars["Fh"] = 2 # f_scaled (for v.1.1)
pars["ENV"] = 500 # Units: L
pars["r"] = 0.25   # Units: day-1
pars["step"] = 1  # Units: day
pars["epsilon"] = 20 # Units: L host-1, day-1 (Rounded estimate from Civitello and Rohr)
pars["sigma"] = 0.5 
pars["m_M"] = 1   # Units: day-1
pars["m_Z"] = 1   # Units: day-1
pars["M_in"] = 10 # Units: day-1
pars["K"] = 5
pars["Det"] = 0 # Units mg C L-1 d-1 (detritus)
pars["DV"] = 0.0096 # Units: mg mm^-3 (Empirical weight - length regression; note not units of carbon)
# This correction factor is explained in supplement of Civitello et al. 2020 Proc B, based on repeated 2 vs 22 hour sheds by RBH
pars["yRP"] = pars["yRP"]*17.5 # new yRP = 0.824 * 17.5 (1-8-19)( # old yRP = 0.0471 * (1/0.4) * 7 [expected shedding output per weekly shed * (snails shed 40% of their total cercariae during 9-11 AM.) * 7 days] 
pars["hb"] = 0.001
pars["Hatch"] = 0.5

######

# DEB() simulations one time step of a DEB model for one snail in a potentially shared environment
DEB = function(step, Food, L, e, D, RH, P, RP, DAM, HAZ, iM, k, M, EM, 
               Fh, muD, DR, yRP, ph, yPE, iPM, eh, mP, alpha, yEF,
               LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp, SAtotal, r, K, Det){
  # starting conditions 
  initials = c(Food=Food, L=L, e=e, D=D, RH=RH, P=P, RP=RP, DAM=DAM, HAZ=HAZ)
  # deb parameters
  parameters = c(iM, k, M, EM, Fh, muD, DR, yRP, ph, yPE, iPM,
                 eh, mP, alpha, yEF, LM, kd, z, kk, hb, theta, mR, yVE, ENV, Lp, SAtotal, r, K, Det)
  # estimate starting deb conditions using fitted params by solving ode's
  ## return survival and host shell length  
  DEBstep <- lsoda(initials, c(0,step), func = "derivs", dllname = "IndividualModel_IBM3", 
                   initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), maxsteps=500000,
                   as.numeric(parameters),  rtol=1e-6, atol=1e-6, hmax=1)
  DEBstep[2, 2:12] # 12 = survival
}

Initialize_IBM = function(N.snail, min.L, max.L, Food){
  L = runif(N.snail, min = min.L, max = max.L)
  e = rep(0.9, times=N.snail)
  snail.stats = cbind("ID" = 1:N.snail, "L" = L, "e" = e, "D" = rep(0, times=N.snail), "RH" = rep(0, times=N.snail),
                      "P" = rep(0, times=N.snail), "RP" = rep(0, times=N.snail), "DAM" = rep(0, times=N.snail),
                      "HAZ" = rep(0, times=N.snail), "LG" = L, "DEBmass" = pars["M"]*((1 + e*pars["EM"])/(1 + pars["EM"]))*L^3, 
                      "Appmass" = pars["DV"]*L^3, "repro" = rep(0, times=N.snail), "Cercs" = rep(0, times=N.snail), "t" = rep(0, times=N.snail))
  list("Snails" = list(snail.stats), "Env_F" = Food, "Env_M" = 0, "Env_Z" = 0, "Env_G" = 0)
}



### Exposure submodel
Infection = function(snail.stats, miracidia, parameters){
  # Parameters
  epsilon = as.numeric(parameters["epsilon"])
  sigma = as.numeric(parameters["sigma"])
  ENV = as.numeric(parameters["ENV"])
  m_M = as.numeric(parameters["m_M"])
  step = as.numeric(parameters["step"])
  
  # Later calculations depend on exposure probabilities
  exp.rates =epsilon/ENV*(snail.stats[,"L"]>0) # This is just to get uniform exposure rates
  sum.exp.rates = sum(exp.rates)
  
  # Probabilities for fate of miracidia
  P.left.in.water = exp(-(m_M+sum(exp.rates))*step)                             # Still in water
  P.infects.this.snail = (1 - P.left.in.water)*(sigma*exp.rates/(m_M+sum.exp.rates))  # Infect a snail
  # Die in water or fail to infect
  P.dead = (1 - P.left.in.water)*(m_M/(m_M+sum.exp.rates)) + sum((1 - P.left.in.water)*((1-sigma)*exp.rates/(m_M+sum.exp.rates)))                      # die
  
  prob.vector = c(P.infects.this.snail, P.left.in.water, P.dead)
  
  # Multinomial outcome
  rmultinom(n=1, size=miracidia, prob=prob.vector)
  #sum(P.left.in.water, P.invades.this.snail, P.dead)
}


# Start controlling NetLogo
state = Initialize_IBM(60, 2, 16, 1)
n.ticks = 112
Env_G = numeric()
snail.IDs = length(state$Snails[[1]][,"ID"])

for(t in 1:n.ticks){
  # Daily check for time-varying parameters
  #pars["Det"] = ifelse(t %% pars["pulse_freq"] == 0, 1, 0) 
  pars["M_in"] = ifelse( t %in% c(7, 21, 42, 63), 150, 0)
  
  # Get environmental stats
  environment = c("F" = state$Env_F[t], "M" = state$Env_M[t], "Z" = state$Env_Z[t], "G" = state$Env_G[t])
  
  # Simulate snails
  if(dim(state$Snails[[t]])[1] > 0){
    # set host variables
    snail.stats = state$Snails[[t]]
    N.snails = length(snail.stats[,"L"])

    # Infect snails
    Infection.step = as.vector(Infection(snail.stats, environment[2], pars)) # Who gets infected
    snail.stats[which(Infection.step[1:N.snails] > 0),"P"] = snail.stats[which(Infection.step[1:N.snails] > 0),"P"] + 2.85e-5 # add biomass of one miracidia

    # Update DEBS, HAZ=0 so survival probs are calculated for the current day
    snail.update = t(mapply(DEB, L=snail.stats[,"L"], e=snail.stats[,"e"], D=snail.stats[,"D"], RH=snail.stats[,"RH"],
                            P=snail.stats[,"P"], RP=snail.stats[,"RP"], DAM=snail.stats[,"DAM"], Lp=snail.stats[,"LG"], 
                            MoreArgs = list(step=1, HAZ=0, Food=as.numeric(environment["F"]),
                                            iM=pars["iM"], k=pars["k"], M=pars["M"], EM=pars["EM"], Fh=pars["Fh"], 
                                            muD=pars["muD"],
                                            DR=pars["DR"], yRP=pars["yRP"], ph=pars["ph"], yPE=pars["yPE"], iPM=pars["iPM"], eh=pars["eh"],
                                            mP=pars["mP"], alpha=pars["alpha"], yEF=pars["yEF"], LM=pars["LM"], kd=pars["kd"], z=pars["z"], 
                                            kk=pars["kk"], hb=pars["hb"],theta=pars["theta"], mR=pars["mR"], yVE=pars["yVE"], 
                                            SAtotal= sum(snail.stats[,"L"]^2), 
                                            ENV=pars["ENV"], r=pars["r"], K=pars["K"], 
                                            Det=pars["Det"]))) # detritus (Det) defined in C file
    
    snail.update[is.nan(snail.update)] <- 0 # turn nans from matrix into 0
    
    # Calculate some useful individual-level stuff
    DEBmass = pars["M"]*((1 + snail.update[,"e"]*pars["EM"])/(1 + pars["EM"]))*snail.update[,"L"]^3
    Appmass = pars["DV"]*snail.update[,"LG"]^3
    
    Eggs = floor(snail.update[,"RH"]/0.015)  # Figure out how many (whole) eggs are released
    snail.update[,"RH"] = snail.update[,"RH"] %% 0.015        # Remove released cercariae from the buffer
    Cercs = floor(snail.update[,"RP"]/4e-5)  # Figure out how many (whole) cercs are released
    snail.update[,"RP"] = snail.update[,"RP"] %% 4e-5         # Remove released cercariae from buffer
    
    # Store calculated individual-level stuff
    snail.update = cbind(snail.update, "DEBmass" = DEBmass, "Appmass" = Appmass, "repro" = Eggs, "Cercs" = Cercs, "t" = rep(t, times=length(Cercs)))
    
    # Snails have to survive the day
    mortality.luck = runif(N.snails, min=0, max=1)
    surviving = mortality.luck >= 1 - exp(-snail.update[,"HAZ"])
    state$Snails[[t+1]] = cbind("ID" = snail.stats[which(surviving == TRUE), "ID"],  snail.update[which(surviving == TRUE),c("L", "e", "D", "RH", "P", "RP", "DAM", "HAZ", "LG", "DEBmass", "Appmass", "repro", "Cercs", "t")])
    
    # Update environment
    state$Env_F[t+1] = max(0, snail.update[1,"Food"])
    state$Env_M[t+1] = as.numeric(Infection.step[N.snails + 1] + pars["M_in"]) # total miracidia density 
    state$Env_Z[t+1] = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"]) + sum(Cercs)/pars["ENV"]) # total cerc density

    }else{ # if there are no hosts  
    
    state$Env_M[t+1] = as.numeric(environment[2]*exp(-pars["m_M"]*pars["step"]) + pars["M_in"]) # total miracidia density 
    state$Env_Z[t+1] = as.numeric(environment[3]*exp(-pars["m_Z"]*pars["step"])) # total cerc density
    state$Env_F[t+1] = ifelse(pars["Det"] == 0, as.numeric(pars["K"]*environment[1]/(environment[1] + (pars["K"] - environment[1])*exp(-pars["r"]*pars["step"]))), as.numeric(environment[1] + pars["Det"]))
    Eggs <- 0
    }
    
  state$Env_G[t+1] = max(0, sum(Eggs))

  if(t > 10){
    births =  rbinom(n=1, size=state$Env_G[t - 10], prob=pars["Hatch"])
    if(births > 0){
      neonates = cbind("ID" = (snail.IDs+1):(snail.IDs + births), "L" = rep(0.75, births), "e" = rep(0.9, times=births), "D" = rep(0, times=births), "RH" = rep(0, times=births),
                  "P" = rep(0, times=births), "RP" = rep(0, times=births), "DAM" = rep(0, times=births),
                  "HAZ" = rep(0, times=births), "LG" = rep(0.75, births), "DEBmass" = pars["M"]*((1 + 0.9*pars["EM"])/(1 + pars["EM"]))*0.75^3,
                  "Appmass" = pars["DV"]*0.75^3, "repro" = rep(0, times=births), "Cercs" = rep(0, times=births), "t" = rep(t, times=births))
      snail.IDs = snail.IDs + births # Updates the running count of snails
      state$Snails[[t+1]] = rbind(state$Snails[[t+1]], neonates)
    }
  }
  print(t)
}

digest_state = function(state){
  span = length(state$Env_F)
  N.snails = unlist(lapply(state$Snails, FUN=length))/dim(state$Snails[[1]])[2]
  DEBmass = matrix(unlist(lapply(state$Snails, FUN=colSums)), ncol=dim(state$Snails[[1]])[2], byrow=T)[, which(colnames(state$Snails[[1]]) == "DEBmass")]
  Appmass = matrix(unlist(lapply(state$Snails, FUN=colSums)), ncol=dim(state$Snails[[1]])[2], byrow=T)[, which(colnames(state$Snails[[1]]) == "Appmass")]
  data.frame("time" = 1:span, "Snails" = N.snails, "DEBmass" = DEBmass, "Appmass" = Appmass, "Food" = state$Env_F,
             "Miracidia" = state$Env_M, "Cercariae" = state$Env_Z, "Eggs" = state$Env_G)
}



plot(digest_state(state))

digest_i_states = function(state){
  all_inds = do.call(rbind, state$Snails)
  inf_inds = subset(data.frame(all_inds), Cercs > 0)
  inf_inds
}

pop.stats = digest_state(state)
inf.snails = digest_i_states(state)

library(ggplot2)
p.cercs = ggplot(data=inf.snails, aes(x=t, y=Cercs, group=ID, color=ID)) + geom_line()
p.Snails = ggplot(data=pop.stats, aes(x=time, y=Snails)) + geom_line()

snail_avgs = matrix(unlist(lapply(state$Snails, FUN=colMeans)), nrow=length(state$Snails), ncol=dim(state$Snails[[1]])[2], byrow=T)
colnames(snail_avgs) = c("ID", "L", "e", "D", "RH", "P", "RP", "DAM", "HAZ", "LG", "DEBmass", "Appmass", "repro", "Cercs", "t")
snail_avgs
plot(0:n.ticks, snail_avgs[,"e"], typ="l")
