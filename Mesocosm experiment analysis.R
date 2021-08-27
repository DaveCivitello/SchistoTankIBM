# Analysis of tank dynamics and individual-level shedding data from the 2016 mesocosm experiment

### Variables and factor-smooth interactions are set up according to
### advice from http://www.sfs.uni-tuebingen.de/~jvanrij/Tutorial/GAMM.html
### advice from conversation with Lance Waller, RE: using unordered factors and proper variable coding


# 0 - individual-level raw data
setwd("C:/RData")
individuals = read.csv("2016_Tank_master.csv")
founders = subset(individuals, Week == 0)
hist(founders$Diameter)

small_tanks = c(4, 5, 6, 9, 10, 14, 22, 25, 26, 28, 30, 32)
intermediate_tanks = c(2, 3, 13, 21, 24, 33, 34, 37, 40, 41, 43, 46)
large_tanks = c(8, 11, 17, 18, 19, 27, 31, 36, 39, 42, 44, 45)
inf_tanks = c(small_tanks, intermediate_tanks, large_tanks)

inf_founders = subset(founders, Tank %in% inf_tanks)
inf_founders[,"Size"] = ifelse(inf_founders[,"Tank"] %in% small_tanks, "Small", 
                               ifelse(inf_founders[,"Tank"] %in%  intermediate_tanks, "Intermediate", "Large"))
founder_means = aggregate(Diameter ~ Size, FUN=mean, data=inf_founders)

SEM = function(x){
  sd(x)/sqrt(na.omit(length(x)))
}

founder_SDs = aggregate(Diameter ~ Size, FUN=sd, data=inf_founders)

founder_summary = data.frame(founder_means, "SD" = founder_SDs$Diameter)
founder_summary

# 1 - tank-level raw data
tanks = read.csv("Schisto_tanks_2016_for_GAMMs.csv")
snails = read.csv("Schisto_tanks_2016_Individual_Infections.csv")
tanks[,"Tank"] = as.factor(tanks[,"Tank"]) # Makes sure tank is ready for random effect
tanks = subset(tanks, Schisto=="Yes" & Week != 0) # Focus on 36 tanks that received schisto
tanks[,"High"] = as.numeric(tanks[,"N...P"] == "High") # Recode as "high fertilizer" 1 or 0
tanks[,"Size"] = factor(tanks[,"Size"], levels=c("Small", "Intermediate", "Large"), ordered=F)
tanks[,"PercapShed"] = as.numeric(tanks[,"Cercarial_production"] / tanks[,"Infected_abundance"])

# 2 - treatment level summary statistics for plotting



# Snail density
density_m = aggregate(Snail_density ~ High*Size*Week, data=tanks, FUN=mean)
density_se = aggregate(Snail_density ~ High*Size*Week, data=tanks, FUN=SEM)
# Snail biomass
biomass_m = aggregate(Biomass_density ~ High*Size*Week, data=tanks, FUN=mean)
biomass_se = aggregate(Biomass_density ~ High*Size*Week, data=tanks, FUN=SEM)
# Periphyton
peri_m = aggregate(Peri_F ~ High*Size*Week, data=tanks, FUN=mean)
peri_se = aggregate(Peri_F ~ High*Size*Week, data=tanks, FUN=SEM)
# Infected snails
Inf_m = aggregate(Infected_abundance ~ High*Size*Week, data=tanks, FUN=mean, drop=F)
Inf_se = aggregate(Infected_abundance ~ High*Size*Week, data=tanks, FUN=SEM, drop=F)
# Cercariae
cerc_m = aggregate(Cercarial_production ~ High*Size*Week, data=tanks, FUN=mean, drop=F)
cerc_se = aggregate(Cercarial_production ~ High*Size*Week, data=tanks, FUN=SEM, drop=F)
# Per cap shedding
percap_m = aggregate(PercapShed ~ High*Size*Week, data=tanks, FUN=mean, drop=F)
percap_se = aggregate(PercapShed ~ High*Size*Week, data=tanks, FUN=SEM, drop=F)

tanksummary = data.frame(density_m, "Snail_density_SE" = density_se$Snail_density,
                         "Biomass_density" = biomass_m$Biomass_density, "Biomass_density_SE" = biomass_se$Biomass_density,
                         "Periphyton" = peri_m$Peri_F, "Periphyton_SE" = peri_se$Peri_F,
                         "Infected" = Inf_m$Infected_abundance, "Infected_SE" = Inf_se$Infected_abundance,
                         "Cercariae" = cerc_m$Cercarial_production, "Cercariae_SE" = cerc_se$Cercarial_production,
                         "PercapShed" = c(rep(NA, times=18), percap_m$PercapShed, rep(NA, times=12)), "PercapShed_SE" = c(rep(NA, times=18), percap_se$PercapShed, rep(NA, times=12)))

tanksummary[,"High"] = as.factor(tanksummary[,"High"])

library(mgcv)
library(glmmTMB)
library(ggeffects)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(gridExtra)
library(itsadug)
library(gratia)

# 3 - Analysis

### Setting the plotting theme that I like
theme_update(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.position="None", axis.ticks.length = unit(-1.5, "mm"), axis.text.x = element_text(margin = margin(t=6)),
             axis.text.y = element_text(margin = margin(r=6)), axis.title= element_blank())


### Use GAMMs with autocorrelated error to analyze population level dynamics

## Figure 2 - What are the temporal dynamics of snails, snail biomass, and periphyton?

# 2A, B - Snail density
m_aS = gamm(Snail_density ~ Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"), 
            family=quasipoisson, correlation=corCAR1(form=~Week|Tank), data=tanks)
summary(m_aS$gam)
draw(derivatives(m_aS$gam, term = "Week", interval="simultaneous")) # Peaks occur when first derivatives pass from positive to negative
fd_aS = derivatives(m_aS$gam, term = "s(Week):SizeSmall", interval="simultaneous", n=1400)
fd_aS[,"sig_fd"] = ifelse(fd_aS[,"lower"]*fd_aS[,"upper"] > 0, 1, 0)
fd_aS
p2ABfit = plot_smooth(m_aS$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, se=1, shade=T, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p2ABfit[,"High"] = as.factor(p2ABfit[,"High"])

p2A =  ggplot(data=p2ABfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue"), na.value=NA) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red"), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Snail_density, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Snail_density , ymin=Snail_density - Snail_density_SE,
                                       ymax=Snail_density + Snail_density_SE, colour=interaction(Size, High)), inherit.aes = F) 


p2B =  ggplot(data=p2ABfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA), na.value=NA) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Snail_density , colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Snail_density , ymin=Snail_density - Snail_density_SE,
                                       ymax=Snail_density + Snail_density_SE, colour=interaction(Size, High)), inherit.aes = F)

# 2C, D Snail biomass 
m_aS2 = gamm(pmax(Biomass_density,1) ~ Size  + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"), 
             family=Gamma(link="log"), correlation=corCAR1(form=~Week|Tank), niterPQL=100, method="REML", data=tanks)
summary(m_aS2$gam)
#draw(derivatives(m_aS2$gam, term = "Week")) # Peaks occur when first derivatives pass from positive to negative

p2CDfit = plot_smooth(m_aS2$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, se=1, shade=T, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p2CDfit[,"High"] = as.factor(p2CDfit[,"High"])

p2C =  ggplot(data=p2CDfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue"), na.value=NA) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red"), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Biomass_density,colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Biomass_density , ymin=Biomass_density - Biomass_density_SE,
                                       ymax=Biomass_density + Biomass_density_SE, colour=interaction(Size, High)), inherit.aes = F)

p2D =  ggplot(data=p2CDfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA), na.value=NA) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Biomass_density,colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Biomass_density , ymin=Biomass_density - Biomass_density_SE,
                                       ymax=Biomass_density + Biomass_density_SE, colour=interaction(Size, High)), inherit.aes = F)

#2E, F - Periphyton
m_peri = gamm(Peri_F ~  Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"),
              family=Gamma(link="log"), correlation=corCAR1(form=~Week|Tank), data=tanks)
summary(m_peri$gam)

p2EFfit = plot_smooth(m_peri$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, se=1, shade=T,
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p2EFfit[,"High"] = as.factor(p2EFfit[,"High"])

p2E =  ggplot(data=p2EFfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue"), na.value=NA) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red"), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Periphyton, group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Periphyton, ymin=Periphyton - Periphyton_SE,
                                       ymax=Periphyton + Periphyton_SE, colour=interaction(Size, High)), inherit.aes = F)

p2F =  ggplot(data=p2EFfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA), na.value=NA) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Periphyton, group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Periphyton, ymin=Periphyton - Periphyton_SE,
                                       ymax=Periphyton + Periphyton_SE, colour=interaction(Size, High)), inherit.aes = F)


### Create a blank "spacer" for cleaner plotting ###
spacer = ggplot(data=p2ABfit, aes(x=Week, y=fit)) +
  geom_blank() + theme_void()


### Fig 2
Fig2 = plot_grid(spacer, spacer, spacer, spacer,
                 spacer, p2A, p2B, spacer,
                 spacer, p2C, p2D,spacer,
                 spacer, p2E, p2F, spacer,
                 spacer, spacer, spacer, spacer,
                 align="hv", ncol=4, nrow=5, rel_widths=c(0.1, 1, 1, 0.05), rel_heights=c(0.15, 1, 1, 1, 0.1), axis="rltb", scale=1) +
  # Column labels
  draw_label(expression(underline("High nutrient")), x=0.33, y=0.99) +
  draw_label(expression(underline("Low nutrient")), x=0.78, y=0.99) +
  # y-axis labels
  draw_label("Host density ± SE", x=0.02, y=0.81, angle=90, size=14) +
  draw_label("Host biomass density, mg ± SE", x=0.02, y=0.5, angle=90, size=14) +
  draw_label("Periphyton productivity ± SE", x=0.02, y=0.19, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail and schistosome introduction", x = 0.55, y = 0.01, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.12, y = 0.95) +  draw_label("B", x = 0.59, y = 0.95) +
  draw_label("C", x = 0.12, y = 0.64) +  draw_label("D", x = 0.59, y = 0.64) +
  draw_label("E", x = 0.12, y = 0.33) +  draw_label("F", x = 0.59, y = 0.33) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.68")), x=0.9, y=0.94, size=12)+
  draw_label(expression(paste(R^2, " = 0.48")), x=0.9, y=0.62, size=12)+
  draw_label(expression(paste(R^2, " = 0.27")), x=0.9, y=0.32, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.42, y=0.95, size=11) +
  draw_label("Small", x = 0.42, y=0.93, size=11, color = "blue") +
  draw_label("Medium", x = 0.42, y=0.91, size=11, color = "black") +
  draw_label("Large", x = 0.42, y=0.89, size=11, color = "red") 

save_plot("Fig2_SufferComp.png", Fig2, ncol=4, nrow=5, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")



## Figure 3 - infected snails, and cercariae
# 1B 

# Needed to add a single infected snail to one tank from each replicate at week 2 to avoid convergence issue from log-link
tanks[,"Infected_abundance2"] = tanks[,"Infected_abundance"]
tanks[c(37, 39, 41, 42, 45, 51),"Infected_abundance2"] = 1

m_Inf = gamm(Infected_abundance2 ~ Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"), niterPQL=100,
             family=quasipoisson(link="log"), correlation=corCAR1(form=~Week|Tank), data=subset(tanks, Week >= 1))
summary(m_Inf$gam)

p3ABfit = plot_smooth(m_Inf$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", shade=T, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p3ABfit[,"High"] = as.factor(p3ABfit[,"High"])

p3A =  ggplot(data=p3ABfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue"), na.value=NA) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red"), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Infected, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Infected, ymin=Infected - Infected_SE, ymax=Infected + Infected_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)

p3A


p3B = ggplot(data=p3ABfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA), na.value=NA) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Infected, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Infected, ymin=Infected - Infected_SE, ymax=Infected + Infected_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)
p3B

#1C - cercariae
m_Cerc = gamm(Cercarial_production ~ Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"),
              family=quasipoisson, correlation=corCAR1(form=~Week|Tank), data=tanks)
summary(m_Cerc$gam)

draw(derivatives(m_Cerc$gam, term = "Week")) # Peaks occur when first derivatives pass from positive to negative

p3CDfit = plot_smooth(m_Cerc$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, se=1, shade=T,
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p3CDfit[,"High"] = as.factor(p3CDfit[,"High"])

p3C =  ggplot(data=p3CDfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue"), na.value=NA) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red"), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Cercariae, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Cercariae, ymin=Cercariae - Cercariae_SE, ymax=Cercariae + Cercariae_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)

p3D =  ggplot(data=p3CDfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA), na.value=NA) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Cercariae, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Cercariae, ymin=Cercariae - Cercariae_SE, ymax=Cercariae + Cercariae_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)


#### Figure 3, population-level infection/cerc data, splitting things out by nutrient level ##

Fig3 = plot_grid(spacer, spacer, spacer, spacer,
                    spacer, p3A, p3B, spacer,
                    spacer, p3C, p3D,spacer,
                    spacer, spacer, spacer, spacer,
                    align="hv", ncol=4, nrow=4, rel_widths=c(0.1, 1, 1, 0.05), rel_heights=c(0.15, 1, 1, 0.1), axis="rltb", scale=1) +
  # Column labels
  draw_label(expression(underline("High nutrient")), x=0.33, y=0.99) +
  draw_label(expression(underline("Low nutrient")), x=0.78, y=0.99) +
  # y-axis labels
  draw_label("Infected host density ± SE", x=0.02, y=0.75, angle=90, size=14) +
  draw_label("Total cercarial production ± SE", x=0.02, y=0.25, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail and schistosome introduction", x = 0.55, y = 0.01, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.11, y = 0.92) +  draw_label("B", x = 0.57, y = 0.92) +
  draw_label("C", x = 0.11, y = 0.47) +  draw_label("D", x = 0.57, y = 0.47) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.33")), x=0.9, y=0.94, size=12)+
  draw_label(expression(paste(R^2, " = 0.71")), x=0.9, y=0.47, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.42, y=0.95, size=11) +
  draw_label("Small", x = 0.42, y=0.93, size=11, color = "blue") +
  draw_label("Medium", x = 0.42, y=0.91, size=11, color = "black") +
  draw_label("Large", x = 0.42, y=0.89, size=11, color = "red") 

save_plot("Fig3_SufferComp.png", Fig3, ncol=4, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")


## Figure 4 - Population- and Individual-level analyses of per capita shedding?

#4A, B - per capita cercariae, pop level
m_CercPC = gamm(Cercarial_production ~ Size + offset(log(Infected_abundance)) + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"),
                family=quasipoisson, correlation=corCAR1(form=~Week|Tank), data=subset(tanks, Week>=4 & Infected_abundance > 0))
summary(m_CercPC$gam)
draw(derivatives(m_CercPC$gam, term = "Week")) # Peaks occur when first derivatives pass from positive to negative (no peak, always decreasing)

p4ABfit = plot_smooth(m_CercPC$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, se=1, shade=T,
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p4ABfit[,"High"] = as.factor(p4ABfit[,"High"])


p4A =  ggplot(data=p4ABfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue"), na.value=NA) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red"), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=subset(tanksummary, Week >= 4), aes(x=Week, y=PercapShed, group=interaction(Size, High))) +
  geom_linerange(data=subset(tanksummary, Week >= 4), aes(x=Week, y=PercapShed, ymin=PercapShed - PercapShed_SE,
                                                          ymax=PercapShed + PercapShed_SE, colour=interaction(Size, High)), inherit.aes = F)

p4B =  ggplot(data=p4ABfit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA), na.value=NA) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA), na.value=NA) + xlim(c(1, 15)) +
  geom_point(data=subset(tanksummary, Week >= 4), aes(x=Week, y=PercapShed, group=interaction(Size, High))) +
  geom_linerange(data=subset(tanksummary, Week >= 4), aes(x=Week, y=PercapShed, ymin=PercapShed - PercapShed_SE,
                                                          ymax=PercapShed + PercapShed_SE, colour=interaction(Size, High)), inherit.aes = F)
# 4C, D per capita shedding, ind level

# Center Diameter and total biomass to help glmmTMB
snails[,"Comp_Biomass"] = snails[,"Total_Biomass"] - snails[,"Biomass"]
high_tanks = unique(subset(tanks, High == 1, select="Tank"))[,1]

snails[,"High"] = ifelse(snails[,"Tank"] %in% high_tanks, 1, 0)
snails[,"Alg_prod"] = 0
for(i in 1:length(snails[,"Tank"])){
  snails[i,"Alg_prod"] = subset(tanks, Week == snails[i, "Week"] & Tank == snails[i, "Tank"], select="Peri_F")
}

snails2 = transform(snails, Diameter=scale(Diameter,center=T, scale=F), Biomass=scale(Biomass/1000, center=T, scale=F), Total_Biomass=scale(Total_Biomass/1000, center=T, scale=F),
                    Comp_Biomass=scale(Comp_Biomass/1000, center=T, scale=F), Alg_prod = scale(Alg_prod/1000, center=T, scale=F))
mod2 = glmmTMB(Cercariae ~ Biomass*Comp_Biomass + (1|Tank), family=nbinom2, data=snails2)
summary(mod2)

sizes = quantile(snails2$Biomass, probs=c(0.1, 0.9))

small = ggpredict(mod2, terms="Comp_Biomass [-0.53:0.83, by = 0.01]", condition = c(Biomass = as.numeric(sizes[1])))
large = ggpredict(mod2, terms="Comp_Biomass [-0.53:0.83, by = 0.01]", condition = c(Biomass = as.numeric(sizes[2])))

snails2[,"High"] = as.factor(snails2[,"High"])

p4C = ggplot(data=snails2, aes(x=1000*(Comp_Biomass+ 0.52), y=Cercariae)) +
  theme(axis.ticks.length = unit(-1.5, "mm")) + 
  geom_point() + scale_y_log10(limits=c(1, 10000)) +
  geom_line(data=small, aes(1000*(x+0.52), predicted), inherit.aes = F) +
  geom_ribbon(data=small, aes(1000*(x+0.52), predicted, ymin = conf.low, ymax = conf.high), alpha = .2, inherit.aes = F) +
  geom_line(data=large, aes(1000*(x+0.52), predicted), size=2, inherit.aes = F) +
  geom_ribbon(data=large, aes(1000*(x+0.52), predicted, ymin = conf.low, ymax = conf.high), alpha = .2, inherit.aes = F)

p4C

Fig4 = plot_grid(spacer,  spacer, spacer,
                  spacer, p4A, spacer,
                  spacer, p4B,spacer,
                 spacer,  spacer, spacer,
                 spacer, p4C, spacer,
                  spacer, spacer, spacer,
                  align="hv", ncol=3, nrow=6, rel_widths=c(0.075, 1, 0.05), 
                  rel_heights=c(0.05, 1, 1, 0.1, 1, 0.1), axis="rltb", scale=1) +
  # y-axis labels
  draw_label("Per capita cercarial production ± SE", x=0.03, y=0.5, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail introduction", x = 0.55, y = 0.35, size=14) +
  draw_label("Snail biomass density, mg", x = 0.55, y = 0.01, ,size=14) +
  # panel labels
  draw_label("A) High Nutrient", x = 0.16, y = 0.97, hjust=0) + 
  draw_label("B) Low Nutrient", x = 0.16, y = 0.65, hjust=0) +
  draw_label("C", x = 0.16, y = 0.32, hjust = 0) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.67")), x=0.22, y=0.93, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.72, y=0.95, size=11) +
  draw_label("Small", x = 0.72, y=0.93, size=11, color = "blue") +
  draw_label("Medium", x = 0.72, y=0.91, size=11, color = "black") +
  draw_label("Large", x = 0.72, y=0.89, size=11, color = "red") 

save_plot("Fig4.png", Fig4, ncol=3, nrow=6, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")


