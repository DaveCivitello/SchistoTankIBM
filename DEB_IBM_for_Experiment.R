
# Loading packages
library(Matrix)
library(deSolve)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

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
#pars["P"] = 1 # predator density
# pars["PA"] = 1 # predator attack rate units = ENV day-1 (50L sweeped / day or 8.33 * 6 L tanks (from sokolow etal 2014 acta tropica))
# pars["Ph"] = 0.1 # predator handling time; units = day


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
  exp.rates =epsilon/ENV*(snail.stats[,"L"]>0) # This gives uniform exposure rates for all snails
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


run_DEB_IBM = function(N.snail=60, min.L=2, max.L=16, Food=1, n.ticks=112, drip=T){
  
  # Set up initial conditions, data structures
  state = Initialize_IBM(N.snail,min.L, max.L, Food)
  Env_G = numeric()
  snail.IDs = length(state$Snails[[1]][,"ID"])
  
  # run the IBM
  for(t in 1:n.ticks){
    #Daily check for time-varying parameters
    if(drip == F){pars["M_in"] = ifelse( t %in% c(1, 15, 29, 43), 140, 0)}
    
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
      
      # hb <- hb + (pred_a * pred_p) / (1 + pred_a * pred_h * N.snails)
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
  return(state)
}


state = run_DEB_IBM()

### Some functions to process the output ###
catch_snails = function(mat, LG = 2){
  mat[which(mat[,"LG"] >= LG),]
}

catch_state = function(state){
  caught_snails = lapply(state$Snails, FUN=catch_snails)
  list("Snails" = caught_snails, "Env_F" = state$Env_F, "Env_M" = state$Env_M, "Env_Z" = state$Env_Z, "Env_G" = state$Env_G)
}

caught_state = catch_state(state)

digest_state = function(state){
  span = length(state$Env_F)
  N.snails = unlist(lapply(state$Snails, FUN=length))/dim(state$Snails[[1]])[2]
  DEBmass = matrix(unlist(lapply(state$Snails, FUN=colSums)), ncol=dim(state$Snails[[1]])[2], byrow=T)[, which(colnames(state$Snails[[1]]) == "DEBmass")]
  Appmass = matrix(unlist(lapply(state$Snails, FUN=colSums)), ncol=dim(state$Snails[[1]])[2], byrow=T)[, which(colnames(state$Snails[[1]]) == "Appmass")]
  data.frame("time" = 1:span, "Snails" = N.snails, "DEBmass" = DEBmass, "Appmass" = Appmass, "Food" = state$Env_F,
             "Miracidia" = state$Env_M, "Cercariae" = state$Env_Z, "Eggs" = state$Env_G)
}


# Full pairwise plots of pop-level outcomes
plot(digest_state(caught_state))
plot(digest_state(state))


# Finding all infected individuals
digest_i_states = function(state){
  all_inds = do.call(rbind, state$Snails)
  inf_inds = subset(data.frame(all_inds), Cercs > 0)
  inf_inds
}

pop.stats = digest_state(state)
inf.snails = digest_i_states(state)

# A function suitable for replicate() to handle and return the pop-level summarized output
summarized_DEB_IBM = function(N.snail=60, min.L=2, max.L=16, Food=1, n.ticks=112, drip=T){
  digest_state(run_DEB_IBM(N.snail, min.L, max.L, Food, n.ticks, drip))
}

# Gets the mean and 95% CI from a collection of replicated runs
summarized_IBM_reps = function(data, endpoint){
  unlisted_reps = matrix(unlist(data[endpoint,]), nrow=dim(data)[2], ncol=length(data[[endpoint,1]]), byrow=T)
  rep_mean = apply(X=unlisted_reps, MARGIN = 2, FUN=mean)
  rep_CI = t(apply(X=unlisted_reps, MARGIN = 2, FUN=quantile, probs=c(0.025, 0.975)))
  data.frame("time" = 1:length(data[[endpoint,1]]), "Mean" = rep_mean, "CI_L" = rep_CI[,1], "CI_H" = rep_CI[,2])
}

## Make supplement figure for drip vs. pulse of miracidia
# Runs replicate simulations with the same conditions
rep_states_drip = replicate(n=50, summarized_DEB_IBM(), simplify=T)
rep_states_pulse = replicate(n=50, summarized_DEB_IBM(drip=F), simplify=T)

FigS1A_df = data.frame(rbind(summarized_IBM_reps(rep_states_drip, "Snails"), summarized_IBM_reps(rep_states_pulse, "Snails")), "M_in" = rep(c("Constant", "Biweekly pulse"), each=113))
#### Set up means and CIs for supplemental figure ####
pS1A = ggplot(data=FigS1A_df, aes(x=time, y=Mean/500, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  #xlab("Day") + ylab(expression(paste("Snail density,  ", L^-1, "± 95% CI"))) +
  labs(fill="Miracidia input mode", colour="Miracidia input mode") +  theme(legend.position = c(0.5, 0.25)) + 
  geom_line() + geom_ribbon(aes(ymin=CI_L/500, ymax=CI_H/500, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("red", "black")) +
  scale_color_manual(values=c("red", "black"))

pS1A

FigS1B_df = data.frame(rbind(summarized_IBM_reps(rep_states_drip, "Cercariae"), summarized_IBM_reps(rep_states_pulse, "Cercariae")), "M_in" = rep(c("constant", "biweekly pulse"), each=113))
#### Set up means and CIs for supplemental figure ####
pS1B = ggplot(data=FigS1B_df, aes(x=time, y=Mean, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  #xlab("Day") + ylab(expression(paste("Snail density,  ", L^-1, "± 95% CI"))) +  
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L, ymax=CI_H, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("red", "black")) +
  scale_color_manual(values=c("red", "black"))

pS1B

FigS1C_df = data.frame(rbind(summarized_IBM_reps(rep_states_drip, "Food"), summarized_IBM_reps(rep_states_pulse, "Food")), "M_in" = rep(c("constant", "biweekly pulse"), each=113))
#### Set up means and CIs for supplemental figure ####
pS1C = ggplot(data=FigS1C_df, aes(x=time, y=Mean, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  #xlab("Day") + ylab(expression(paste("Snail density,  ", L^-1, "± 95% CI"))) +  
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L, ymax=CI_H, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("red", "black")) +
  scale_color_manual(values=c("red", "black"))

pS1C

FigS1D_df = data.frame(rbind(summarized_IBM_reps(rep_states_drip, "Miracidia"), summarized_IBM_reps(rep_states_pulse, "Miracidia")), "M_in" = rep(c("constant", "biweekly pulse"), each=113))
#### Set up means and CIs for supplemental figure ####
pS1D = ggplot(data=FigS1D_df, aes(x=time, y=Mean/500, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  #xlab("Day") + ylab(expression(paste("Snail density,  ", L^-1, "± 95% CI"))) +  
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L/500, ymax=CI_H/500, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("red", "black")) +
  scale_color_manual(values=c("red", "black"))

pS1D

spacer = ggplot(data=FigS1A_df, aes(x=time, y=Mean)) +
  geom_blank() + theme_void()

FigS1 = plot_grid(spacer, pS1A, spacer, pS1B, 
                  spacer, pS1C, spacer, pS1D,
                  spacer, spacer, spacer, spacer, 
                  nrow=3, ncol=4, rel_widths = c(0.05, 1, 0.05, 1), rel_heights = c(1, 1, 0.05)) +
  # Panel y-axis labels
  draw_label(expression(paste("Snail density,  ", L^-1, "± 95% CI")), angle=90, x=0.02, y=0.75) +
  draw_label(expression(paste("Cercariae density,  ", L^-1, "± 95% CI")), angle=90, x=0.52, y=0.75) +
  draw_label(expression(paste("Resource density, mg C  ", L^-1, "± 95% CI")), angle=90, x=0.02, y=0.25) +
  draw_label(expression(paste("Miracidia density,  ", L^-1, "± 95% CI")), angle=90, x=0.52, y=0.25) +
  # X-axis label
  draw_label("Day", x=0.52, y=0.025) +
  # Panel labels
  draw_label("A", x=0.07, y=0.975) +
  draw_label("B", x=0.59, y=0.975) +
  draw_label("C", x=0.07, y=0.48) +
  draw_label("D", x=0.58, y=0.48) 
  

save_plot("FigS1_IBM.png", FigS1, ncol=2, nrow=2, base_height=4, base_aspect_ratio = 1.1, dpi=600, units="in")


####


## Make supplement figure 2 for gradients of miracidia input and resource productivity
pars["M_in"] = 1
rep_states_drip1 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["M_in"] = 10
rep_states_drip10 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["M_in"] = 100
rep_states_drip100 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["M_in"] = 1000
rep_states_drip1000 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["M_in"] = 10000
rep_states_drip10000 = replicate(n=50, summarized_DEB_IBM(), simplify=T)

pars["M_in"] = 10
pars["r"] = 0.05
rep_states_r1 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["r"] = 0.125
rep_states_r2 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["r"] = 0.25
rep_states_r3 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["r"] = 0.5
rep_states_r4 = replicate(n=50, summarized_DEB_IBM(), simplify=T)
pars["r"] = 1
rep_states_r5 = replicate(n=50, summarized_DEB_IBM(), simplify=T)



d1 = summarized_IBM_reps(rep_states_drip1, "Snails")
d10 = summarized_IBM_reps(rep_states_drip10, "Snails")
d100 = summarized_IBM_reps(rep_states_drip100, "Snails")
d1000 = summarized_IBM_reps(rep_states_drip1000, "Snails")
d10000 = summarized_IBM_reps(rep_states_drip10000, "Snails")

FigS2A_df = data.frame(rbind(d1, d10, d100, d1000, d10000), "M_in" = rep(c(1, 10, 100, 1000, 10000), each=dim(d1)[1]))
FigS2A_df[, "M_in"] = as.factor(FigS2A_df[, "M_in"])

pS2A = ggplot(data=FigS2A_df, aes(x=time, y=Mean/500, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  labs(fill=expression(paste("Miracidia input rate,  ", d^-1)), colour=expression(paste("Miracidia input rate,  ", d^-1))) +
  theme(legend.position = c(0.1, 0.825)) + 
  geom_line() + geom_ribbon(aes(ymin=CI_L/500, ymax=CI_H/500, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f")) +
  scale_color_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f"))

pS2A


d1 = summarized_IBM_reps(rep_states_drip1, "Cercariae")
d10 = summarized_IBM_reps(rep_states_drip10, "Cercariae")
d100 = summarized_IBM_reps(rep_states_drip100, "Cercariae")
d1000 = summarized_IBM_reps(rep_states_drip1000, "Cercariae")
d10000 = summarized_IBM_reps(rep_states_drip10000, "Cercariae")

FigS2C_df = data.frame(rbind(d1, d10, d100, d1000, d10000), "M_in" = rep(c(1, 10, 100, 1000, 10000), each=dim(d1)[1]))
FigS2C_df [, "M_in"] = as.factor(FigS2C_df [, "M_in"])

pS2C = ggplot(data=FigS2C_df, aes(x=time, y=Mean, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L, ymax=CI_H, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f")) +
  scale_color_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f"))

pS2C

d1 = summarized_IBM_reps(rep_states_drip1, "Food")
d10 = summarized_IBM_reps(rep_states_drip10, "Food")
d100 = summarized_IBM_reps(rep_states_drip100, "Food")
d1000 = summarized_IBM_reps(rep_states_drip1000, "Food")
d10000 = summarized_IBM_reps(rep_states_drip10000, "Food")

FigS2E_df = data.frame(rbind(d1, d10, d100, d1000, d10000), "M_in" = rep(c(1, 10, 100, 1000, 10000), each=dim(d1)[1]))
FigS2E_df [, "M_in"] = as.factor(FigS2E_df [, "M_in"])

pS2E = ggplot(data=FigS2E_df, aes(x=time, y=Mean, group=M_in, colour=M_in)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L, ymax=CI_H, colour=NULL, fill=M_in), alpha=0.2) +
  scale_fill_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f")) +
  scale_color_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f"))

pS2E


d1 = summarized_IBM_reps(rep_states_r1, "Snails")
d10 = summarized_IBM_reps(rep_states_r2, "Snails")
d100 = summarized_IBM_reps(rep_states_r3, "Snails")
d1000 = summarized_IBM_reps(rep_states_r4, "Snails")
d10000 = summarized_IBM_reps(rep_states_r5, "Snails")

FigS2B_df = data.frame(rbind(d1, d10, d100, d1000, d10000), "r" = rep(c(0.05, 0.125, 0.25, 0.5, 1), each=dim(d1)[1]))
FigS2B_df[, "r"] = as.factor(FigS2B_df[, "r"])


pS2B = ggplot(data=FigS2B_df, aes(x=time, y=Mean/500, group=r, colour=r)) +
  xlab(NULL) + ylab(NULL) +
  labs(fill=expression(paste("Resource productivity,  ", d^-1)),
       colour=expression(paste("Resource productivity,  ", d^-1))) +
  theme(legend.position = c(0.1, 0.825)) + 
  geom_line() + geom_ribbon(aes(ymin=CI_L/500, ymax=CI_H/500, colour=NULL, fill=r), alpha=0.2) +
  scale_fill_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f")) +
  scale_color_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f"))

pS2B


d1 = summarized_IBM_reps(rep_states_r1, "Cercariae")
d10 = summarized_IBM_reps(rep_states_r2, "Cercariae")
d100 = summarized_IBM_reps(rep_states_r3, "Cercariae")
d1000 = summarized_IBM_reps(rep_states_r4, "Cercariae")
d10000 = summarized_IBM_reps(rep_states_r5, "Cercariae")

FigS2D_df = data.frame(rbind(d1, d10, d100, d1000, d10000), "r" = rep(c(0.05, 0.125, 0.25, 0.5, 1), each=dim(d1)[1]))
FigS2D_df[, "r"] = as.factor(FigS2D_df[, "r"])

pS2D = ggplot(data=FigS2D_df, aes(x=time, y=Mean, group=r, colour=r)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L, ymax=CI_H, colour=NULL, fill=r), alpha=0.2) +
  scale_fill_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f")) +
  scale_color_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f"))

pS2D

d1 = summarized_IBM_reps(rep_states_r1, "Food")
d10 = summarized_IBM_reps(rep_states_r2, "Food")
d100 = summarized_IBM_reps(rep_states_r3, "Food")
d1000 = summarized_IBM_reps(rep_states_r4, "Food")
d10000 = summarized_IBM_reps(rep_states_r5, "Food")

FigS2F_df = data.frame(rbind(d1, d10, d100, d1000, d10000), "r" = rep(c(0.05, 0.125, 0.25, 0.5, 1), each=dim(d1)[1]))
FigS2F_df[, "r"] = as.factor(FigS2F_df[, "r"])


pS2F = ggplot(data=FigS2F_df, aes(x=time, y=Mean, group=r, colour=r)) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = "None") + 
  geom_line() + geom_ribbon(aes(ymin=CI_L, ymax=CI_H, colour=NULL, fill=r), alpha=0.2) +
  scale_fill_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f")) +
  scale_color_manual(values=c("#2166ac", "black", "#d6604d", "#cb181d", "#67001f"))

pS2F

FigS2 = plot_grid(spacer, pS2A, pS2B, 
                  spacer, pS2C, pS2D,
                  spacer, pS2E, pS2F,
                  spacer, spacer, spacer, 
                  nrow=4, ncol=3, align="hv", rel_widths = c(0.05, 1, 1), rel_heights = c(1, 1, 1, 0.05)) +
  # Panel y-axis labels
  draw_label(expression(paste("Snail density,  ", L^-1, "± 95% CI")), angle=90, x=0.02, y=0.85) +
  draw_label(expression(paste("Cercariae density,  ", L^-1, "± 95% CI")), angle=90, x=0.02, y=0.52) +
  draw_label(expression(paste("Resource density, mg C  ", L^-1, "± 95% CI")), angle=90, x=0.02, y=0.2) +
  # X-axis label
  draw_label("Day", x=0.52, y=0.025) +
  # Panel labels
  draw_label("A", x=0.095, y=0.98) +
  draw_label("B", x=0.58, y=0.98) +
  draw_label("C", x=0.095, y=0.65) +
  draw_label("D", x=0.58, y=0.65) +
  draw_label("E", x=0.095, y=0.32) +
  draw_label("F", x=0.58, y=0.32)

save_plot("FigS2_IBM.png", FigS2, ncol=2, nrow=3, base_height=4, base_aspect_ratio = 1.1, dpi=600, units="in")
