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

## 1 - What are the temporal dynamics of snails, snail biomass, and periphyton?

# 1A, Snail density
m_aS = gamm(Snail_density ~ Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"), 
            family=quasipoisson, correlation=corCAR1(form=~Week|Tank), data=tanks)
summary(m_aS$gam)
draw(derivatives(m_aS$gam, term = "Week", interval="simultaneous")) # Peaks occur when first derivatives pass from positive to negative
fd_aS = derivatives(m_aS$gam, term = "s(Week):SizeSmall", interval="simultaneous", n=1400)
fd_aS[,"sig_fd"] = ifelse(fd_aS[,"lower"]*fd_aS[,"upper"] > 0, 1, 0)
fd_aS
p1fit = plot_smooth(m_aS$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", cex.lab=2, cex.axis=1.5, se=1, shade=T,
                    col=c("black", "black", "black", "limegreen", "limegreen", "limegreen"), lty=c(3, 2, 1, 3, 2, 1), lwd=2, legend_plot_all = F,
                    tcl=0.5, h0=NULL, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p1fit[,"High"] = as.factor(p1fit[,"High"])

p1 =  ggplot(data=p1fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=High, linetype=Size)) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=High, colour=NA), alpha=0.2) +
  scale_linetype_manual(values=c(3, 2, 1)) + scale_color_manual(values=c("black", "limegreen")) +
  scale_fill_manual(values=c("black", "limegreen")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Snail_density , group=interaction(Size, High), shape=Size)) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Snail_density , ymin=Snail_density - Snail_density_SE,
                                  ymax=Snail_density + Snail_density_SE, colour=High), inherit.aes = F)



# 1B, Snail biomass 
m_aS2 = gamm(pmax(Biomass_density,1) ~ Size  + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"), 
             family=Gamma(link="log"), correlation=corCAR1(form=~Week|Tank), niterPQL=100, method="REML", data=tanks)
summary(m_aS2$gam)
draw(derivatives(m_aS2$gam, term = "Week")) # Peaks occur when first derivatives pass from positive to negative

p2fit = plot_smooth(m_aS2$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", cex.lab=2, cex.axis=1.5, se=1, shade=T,
                    col=c("black", "black", "black", "limegreen", "limegreen", "limegreen"), lty=c(3, 2, 1, 3, 2, 1), lwd=2, legend_plot_all = F,
                    tcl=0.5, h0=NULL, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p2fit[,"High"] = as.factor(p2fit[,"High"])


p2 =  ggplot(data=p2fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=High, linetype=Size)) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=High, colour=NA), alpha=0.2) +
  scale_linetype_manual(values=c(3, 2, 1)) + scale_color_manual(values=c("black", "limegreen")) +
  scale_fill_manual(values=c("black", "limegreen")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Biomass_density , group=interaction(Size, High), shape=Size)) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Biomass_density , ymin=Biomass_density - Biomass_density_SE,
                                       ymax=Biomass_density + Biomass_density_SE, colour=High), inherit.aes = F)

p2

#2D - Periphyton
m_peri = gamm(Peri_F ~  Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"),
              family=Gamma(link="log"), correlation=corCAR1(form=~Week|Tank), data=tanks)
summary(m_peri$gam)

p5fit = plot_smooth(m_peri$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", cex.lab=2, cex.axis=1.5, se=1, shade=T,
                    col=c("black", "black", "black", "limegreen", "limegreen", "limegreen"), lty=c(3, 2, 1, 3, 2, 1), lwd=2, legend_plot_all = F,
                    tcl=0.5, h0=NULL, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p5fit[,"High"] = as.factor(p5fit[,"High"])

p5 =  ggplot(data=p5fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=High, linetype=Size)) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=High, colour=NA), alpha=0.2) +
  scale_linetype_manual(values=c(3, 2, 1)) + scale_color_manual(values=c("black", "limegreen")) +
  scale_fill_manual(values=c("black", "limegreen")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Periphyton, group=interaction(Size, High), shape=Size)) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Periphyton, ymin=Periphyton - Periphyton_SE,
                                       ymax=Periphyton + Periphyton_SE, colour=High), inherit.aes = F)

spacer = ggplot(data=p1fit, aes(x=Week, y=fit)) +
  geom_blank() + theme_void()

Fig2 = plot_grid(spacer, spacer, spacer,
                 spacer, p1, spacer,
                 spacer, p2, spacer,
                 spacer, p5, spacer,
                 spacer, spacer, spacer,
                 align="hv", ncol=3, nrow=5, rel_widths=c(0.15, 1, 0.05), rel_heights=c(0.1, 1, 1, 1, 0.1), axis="rltb", scale=1) +
  # y-axis labels
  draw_label("Host density ± SE", x=0.04, y=0.81, angle=90, size=14) +
  draw_label("Host biomass density, mg ± SE", x=0.04, y=0.5, angle=90, size=14) +
  draw_label("Periphyton productivity ± SE", x=0.04, y=0.19, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail and schistosome introduction", x = 0.55, y = 0.01, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.22, y = 0.95) + 
  draw_label("B", x = 0.22, y = 0.64) +
  draw_label("C", x = 0.22, y = 0.33) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.68")), x=0.35, y=0.94, size=12)+
  draw_label(expression(paste(R^2, " = 0.48")), x=0.35, y=0.62, size=12)+
  draw_label(expression(paste(R^2, " = 0.27")), x=0.35, y=0.32, size=12)+
# legend
  draw_label(expression(underline("Founder size")), x=0.86, y=0.95, size=11) +
  draw_label("●", x=0.81, y=0.93, size=16, colour="black") + # unicode character, found using character map
  draw_label("▲", x=0.81, y=0.91, size=12, colour="black") + # unicode character
  draw_label("▪", x=0.81, y=0.89, size=24, colour="black") + # unicode character
  draw_line(x = c(0.825, 0.86), y=c(0.9275, 0.9275), linetype="dotted") + draw_label("Small", x = 0.87, y=0.93, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.9075, 0.9075), linetype="dashed") + draw_label("Medium", x = 0.87, y=0.91, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.8875, 0.8875), linetype="solid") + draw_label("Large", x = 0.87, y=0.89, size=11, hjust=0) +
  draw_label(expression(underline("Nutrients")), x=0.86, y=0.86, size=11) +
  draw_label("High", x=0.86, y=0.84, size=11, colour="limegreen") +
  draw_label("Low", x=0.86, y=0.82, size=11) 

save_plot("Fig2_SufferComp.png", Fig2, ncol=3, nrow=5, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")


### New Fig 2, split out by nutrient treatment
p1H =  ggplot(data=p1fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue")) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Snail_density, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Snail_density , ymin=Snail_density - Snail_density_SE,
                                       ymax=Snail_density + Snail_density_SE, colour=interaction(Size, High)), inherit.aes = F) 


p1L =  ggplot(data=p1fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
   scale_color_manual(values=c("black", NA, "red", NA, "blue", NA)) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA)) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Snail_density , colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Snail_density , ymin=Snail_density - Snail_density_SE,
                                       ymax=Snail_density + Snail_density_SE, colour=interaction(Size, High)), inherit.aes = F)

p2H =  ggplot(data=p2fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue")) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Biomass_density,colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Biomass_density , ymin=Biomass_density - Biomass_density_SE,
                                       ymax=Biomass_density + Biomass_density_SE, colour=interaction(Size, High)), inherit.aes = F)

p2L =  ggplot(data=p2fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA)) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA)) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Biomass_density,colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Biomass_density , ymin=Biomass_density - Biomass_density_SE,
                                       ymax=Biomass_density + Biomass_density_SE, colour=interaction(Size, High)), inherit.aes = F)

p5H =  ggplot(data=p5fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue")) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Periphyton, group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Periphyton, ymin=Periphyton - Periphyton_SE,
                                       ymax=Periphyton + Periphyton_SE, colour=interaction(Size, High)), inherit.aes = F)

p5L =  ggplot(data=p5fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA)) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA)) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Periphyton, group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Periphyton, ymin=Periphyton - Periphyton_SE,
                                       ymax=Periphyton + Periphyton_SE, colour=interaction(Size, High)), inherit.aes = F)

Fig2.v2 = plot_grid(spacer, spacer, spacer, spacer,
                 spacer, p1H, p1L, spacer,
                 spacer, p2H, p2L,spacer,
                 spacer, p5H, p5L, spacer,
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

save_plot("Fig2_SufferComp_v2.png", Fig2.v2, ncol=4, nrow=5, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")



## infected snails, and cercariae?
# 1B 

# Needed to add a single infected snail to one tank from each replicate at week 2 to avoid convergence issue from log-link
tanks[,"Infected_abundance2"] = tanks[,"Infected_abundance"]
tanks[c(37, 39, 41, 42, 45, 51),"Infected_abundance2"] = 1

m_Inf = gamm(Infected_abundance2 ~ Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"), niterPQL=100,
             family=quasipoisson(link="log"), correlation=corCAR1(form=~Week|Tank), data=subset(tanks, Week >= 1))
summary(m_Inf$gam)

p3fit = plot_smooth(m_Inf$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", cex.lab=2, cex.axis=1.5, se=1, shade=T,
                    col=c("black", "black", "black", "limegreen", "limegreen", "limegreen"), lty=c(3, 2, 1, 3, 2, 1), lwd=2, legend_plot_all = F,
                    tcl=0.5, h0=NULL, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p3fit[,"High"] = as.factor(p3fit[,"High"])

p3 =  ggplot(data=p3fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=High, linetype=Size)) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=High, colour=NA), alpha=0.2) +
  scale_linetype_manual(values=c(3, 2, 1)) + scale_color_manual(values=c("black", "limegreen")) +
  scale_fill_manual(values=c("black", "limegreen")) + xlim(c(1, 15)) +
  geom_point(data=subset(tanksummary, Week >= 1), aes(x=Week, y=Infected, group=interaction(Size, High), shape=Size)) +
  geom_linerange(data=subset(tanksummary, Week >= 1), aes(x=Week, y=Infected, ymin=Infected - Infected_SE, ymax=Infected + Infected_SE, colour=High), inherit.aes = F)

p3

#1C - cercariae
m_Cerc = gamm(Cercarial_production ~ Size + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"),
              family=quasipoisson, correlation=corCAR1(form=~Week|Tank), data=tanks)
summary(m_Cerc$gam)

draw(derivatives(m_Cerc$gam, term = "Week")) # Peaks occur when first derivatives pass from positive to negative (no peak, always decreasing)

p4fit = plot_smooth(m_Cerc$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", cex.lab=2, cex.axis=1.5, se=1, shade=T,
                    col=c("black", "black", "black", "limegreen", "limegreen", "limegreen"), lty=c(3, 2, 1, 3, 2, 1), lwd=2, legend_plot_all = F,
                    tcl=0.5, h0=NULL, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p4fit[,"High"] = as.factor(p4fit[,"High"])

p4 =  ggplot(data=p4fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=High, linetype=Size)) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=High, colour=NA), alpha=0.2) +
  scale_linetype_manual(values=c(3, 2, 1)) + scale_color_manual(values=c("black", "limegreen")) +
  scale_fill_manual(values=c("black", "limegreen")) + xlim(c(1, 15)) +
  geom_point(data=subset(tanksummary, Week >= 1), aes(x=Week, y=Cercariae, group=interaction(Size, High), shape=Size)) +
  geom_linerange(data=subset(tanksummary, Week >= 1), aes(x=Week, y=Cercariae, ymin=Cercariae - Cercariae_SE,
                                                          ymax=Cercariae + Cercariae_SE, colour=High), inherit.aes = F)


Fig3A = plot_grid(spacer, spacer, spacer,
                 spacer, p3, spacer,
                 spacer, p4, spacer,
                 spacer, spacer, spacer,
                 align="hv", ncol=3, nrow=4, rel_widths=c(0.15, 1, 0.05), rel_heights=c(0.1, 1, 1, 0.1), axis="rltb", scale=1) +
  # y-axis labels
  draw_label("Infected host density ± SE", x=0.04, y=0.75, angle=90, size=14) +
  draw_label("Total Cercarial production ± SE", x=0.04, y=0.25, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail and schistosome introduction", x = 0.55, y = 0.01, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.20, y = 0.94) + 
  draw_label("B", x = 0.20, y = 0.48) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.30")), x=0.24, y=0.90, size=12)+
  draw_label(expression(paste(R^2, " = 0.67")), x=0.24, y=0.44, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.86, y=0.95, size=11) +
  draw_label("●", x=0.81, y=0.93, size=16, colour="black") + # unicode character, found using character map
  draw_label("▲", x=0.81, y=0.91, size=12, colour="black") + # unicode character
  draw_label("▪", x=0.81, y=0.89, size=24, colour="black") + # unicode character
  draw_line(x = c(0.825, 0.86), y=c(0.9275, 0.9275), linetype="dotted") + draw_label("Small", x = 0.87, y=0.93, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.9075, 0.9075), linetype="dashed") + draw_label("Medium", x = 0.87, y=0.91, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.8875, 0.8875), linetype="solid") + draw_label("Large", x = 0.87, y=0.89, size=11, hjust=0) +
  draw_label(expression(underline("Nutrients")), x=0.86, y=0.86, size=11) +
  draw_label("High", x=0.86, y=0.84, size=11, colour="limegreen") +
  draw_label("Low", x=0.86, y=0.82, size=11) 

save_plot("Fig3A_SufferComp.png", Fig3A, ncol=3, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")



#1C - cercariae
m_CercPC = gamm(Cercarial_production ~ Size + offset(log(Infected_abundance)) + s(Week, by=Size) + s(Week, by=High) + s(Tank, bs="re"),
              family=quasipoisson, correlation=corCAR1(form=~Week|Tank), data=subset(tanks, Week>=4 & Infected_abundance > 0))
summary(m_CercPC$gam)
draw(derivatives(m_CercPC$gam, term = "Week")) # Peaks occur when first derivatives pass from positive to negative (no peak, always decreasing)

p6fit = plot_smooth(m_CercPC$gam, view="Week", plot_all=c("Size", "High"), rug=F, transform=exp, main="", cex.lab=2, cex.axis=1.5, se=1, shade=T,
                    col=c("black", "black", "black", "limegreen", "limegreen", "limegreen"), lty=c(3, 2, 1, 3, 2, 1), lwd=2, legend_plot_all = F,
                    tcl=0.5, h0=NULL, 
                    ylab="Response", rm.ranef=T, xlab="Time (weeks)", hide.label=T)$fv
p6fit[,"High"] = as.factor(p6fit[,"High"])

inf_only = subset(tanks, Week>=4 & Infected_abundance > 0)
plot(PercapShed ~ Week, data=inf_only, pch=21, bg="red")


p6 =  ggplot(data=p6fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=High, linetype=Size)) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=High, colour=NA), alpha=0.2) + 
  scale_linetype_manual(values=c(3, 2, 1)) + scale_color_manual(values=c("black", "limegreen")) +
   scale_fill_manual(values=c("black", "limegreen")) + xlim(c(4, 15)) +
   geom_point(data=subset(tanksummary, Week >= 4), aes(x=Week, y=PercapShed, group=interaction(Size, High), shape=Size)) +
   geom_linerange(data=subset(tanksummary, Week >= 4), aes(x=Week, y=PercapShed, ymin=PercapShed - PercapShed_SE,
                                                           ymax=PercapShed + PercapShed_SE, colour=High), inherit.aes = F)

## 3 - What factors influence per capita shedding?
# Center Diameter and total biomass to help glmmTMB
snails2 = transform(snails, Diameter=scale(Diameter,center=T, scale=F), Total_Biomass=scale(Total_Biomass, center=T, scale=F))
mod2 = glmmTMB(Cercariae ~ Diameter*Total_Biomass + (1|Tank), family=nbinom2, data=snails2)
summary(mod2)
# Body size increases, population (biomass) density decreases shedding

quantile(snails$Diameter, probs=c(0, 0.5, 1))

D3 = ggpredict(mod2, terms="Total_Biomass", condition = c(Diameter = 3 - 13.21405))
D14 = ggpredict(mod2, terms="Total_Biomass", condition = c(Diameter = 14 - 13.21405))
D22 = ggpredict(mod2, terms="Total_Biomass", condition = c(Diameter = 22- 13.21405))

p7 = ggplot(data=snails, aes(x=Total_Biomass, y=Cercariae)) +
  theme(axis.ticks.length = unit(-1.5, "mm")) + 
  geom_point() + scale_y_log10() +
  geom_line(data=D3, aes(x+594.6598, predicted), col="blue") +
  geom_ribbon(data=D3, aes(x+594.6598, predicted, ymin = conf.low, ymax = conf.high), alpha = .2, fill="blue") +
  geom_line(data=D14, aes(x+594.6598, predicted)) +
  geom_ribbon(data=D14, aes(x+594.6598, predicted, ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_line(data=D22, aes(x+594.6598, predicted), col="red") +
  geom_ribbon(data=D22, aes(x+594.6598, predicted, ymin = conf.low, ymax = conf.high), alpha = .2, fill="red")



Fig3B = plot_grid(spacer, spacer, spacer, spacer, spacer,
                  spacer, p3, spacer, p4, spacer,
                  spacer, spacer, spacer, spacer, spacer,
                  spacer, p6, spacer, p7, spacer,
                  spacer, spacer, spacer, spacer, spacer,
                  align="hv", ncol=5, nrow=5, rel_widths=c(0.075, 1, 0.075, 1, 0.05), 
                  rel_heights=c(0.05, 1, 0.1, 1, 0.1), axis="rltb", scale=1) +
  # y-axis labels
  draw_label("Infected host density ± SE", x=0.03, y=0.75, angle=90, size=14) +
  draw_label("Total cercarial production ± SE", x=0.5, y=0.75, angle=90, size=14) +
  draw_label("Per capita cercarial production ± SE", x=0.03, y=0.25, angle=90, size=14) +
  draw_label("Per capita cercarial production ± SE", x=0.5, y=0.25, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail introduction", x = 0.55, y = 0.51, vjust=0 ,size=14) +
  draw_label("Weeks post snail introduction", x = 0.27, y = 0.02, vjust=0 ,size=14) +
  draw_label("Snail biomass density, mg", x = 0.78, y = 0.02, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.12, y = 0.95) + 
  draw_label("B", x = 0.62, y = 0.95) +
  draw_label("C", x = 0.12, y = 0.48) +
  draw_label("D", x = 0.62, y = 0.48) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.30")), x=0.24, y=0.93, size=12)+
  draw_label(expression(paste(R^2, " = 0.66")), x=0.7, y=0.93, size=12)+
  draw_label(expression(paste(R^2, " = 0.67")), x=0.24, y=0.47, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.86, y=0.95, size=11) +
  draw_label("●", x=0.81, y=0.93, size=16, colour="black") + # unicode character, found using character map
  draw_label("▲", x=0.81, y=0.91, size=12, colour="black") + # unicode character
  draw_label("▪", x=0.81, y=0.89, size=24, colour="black") + # unicode character
  draw_line(x = c(0.825, 0.86), y=c(0.9275, 0.9275), linetype="dotted") + draw_label("Small", x = 0.87, y=0.93, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.9075, 0.9075), linetype="dashed") + draw_label("Medium", x = 0.87, y=0.91, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.8875, 0.8875), linetype="solid") + draw_label("Large", x = 0.87, y=0.89, size=11, hjust=0) +
  draw_label(expression(underline("Nutrients")), x=0.86, y=0.86, size=11) +
  draw_label("High", x=0.86, y=0.84, size=11, colour="limegreen") +
  draw_label("Low", x=0.86, y=0.82, size=11) 

save_plot("Fig3B_SufferComp.png", Fig3B, ncol=3, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")

#### Figure 3C, population-level infection/cerc data, splitting things out by nutrient level ##

ggplot(data=p1fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue")) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Snail_density, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Snail_density , ymin=Snail_density - Snail_density_SE,
                                       ymax=Snail_density + Snail_density_SE, colour=interaction(Size, High)), inherit.aes = F) 


p3H =  ggplot(data=p3fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue")) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Infected, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Infected, ymin=Infected - Infected_SE, ymax=Infected + Infected_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)

p3H


p3L = ggplot(data=p3fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA)) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA)) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Infected, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Infected, ymin=Infected - Infected_SE, ymax=Infected + Infected_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)
p3L

p4H =  ggplot(data=p4fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c(NA, "black", NA, "red", NA, "blue")) +
  scale_fill_manual(values=c(NA, NA, NA, "blue", "black", "red")) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Cercariae, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Cercariae, ymin=Cercariae - Cercariae_SE, ymax=Cercariae + Cercariae_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)

p4L =  ggplot(data=p4fit, aes(x=Week, y=fit, group=interaction(Size, High), colour=interaction(Size, High))) +
  geom_line() + geom_ribbon(aes(ymin=ll, ymax=ul, fill=interaction(Size, High), colour=NA), alpha=0.2) +
  scale_color_manual(values=c("black", NA, "red", NA, "blue", NA)) +
  scale_fill_manual(values=c("blue", "black", "red", NA, NA, NA)) + xlim(c(1, 15)) +
  geom_point(data=tanksummary, aes(x=Week, y=Cercariae, colour=interaction(Size, High), group=interaction(Size, High))) +
  geom_linerange(data=tanksummary, aes(x=Week, y=Cercariae, ymin=Cercariae - Cercariae_SE, ymax=Cercariae + Cercariae_SE,
                                       colour=interaction(Size, High)), inherit.aes = F)

Fig3.v2 = plot_grid(spacer, spacer, spacer, spacer,
                    spacer, p3H, p3L, spacer,
                    spacer, p4H, p4L,spacer,
                    spacer, spacer, spacer, spacer,
                    align="hv", ncol=4, nrow=4, rel_widths=c(0.1, 1, 1, 0.05), rel_heights=c(0.15, 1, 1, 0.1), axis="rltb", scale=1) +
  # Column labels
  draw_label(expression(underline("High nutrient")), x=0.33, y=0.99) +
  draw_label(expression(underline("Low nutrient")), x=0.78, y=0.99) +
  # y-axis labels
  draw_label("Infected host density ± SE", x=0.02, y=0.75, angle=90, size=14) +
  draw_label("Total cercarial production ± SE", x=0.02, y=0.35, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail and schistosome introduction", x = 0.55, y = 0.01, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.13, y = 0.92) +  draw_label("B", x = 0.59, y = 0.92) +
  draw_label("C", x = 0.13, y = 0.47) +  draw_label("D", x = 0.59, y = 0.47) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.33")), x=0.9, y=0.94, size=12)+
  draw_label(expression(paste(R^2, " = 0.71")), x=0.9, y=0.47, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.42, y=0.95, size=11) +
  draw_label("Small", x = 0.42, y=0.93, size=11, color = "blue") +
  draw_label("Medium", x = 0.42, y=0.91, size=11, color = "black") +
  draw_label("Large", x = 0.42, y=0.89, size=11, color = "red") 

save_plot("Fig3_SufferComp_v2.png", Fig3.v2, ncol=4, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")


Fig4 = plot_grid(spacer, spacer, spacer, spacer, spacer,
                  spacer, p3, spacer, p4, spacer,
                  spacer, spacer, spacer, spacer, spacer,
                  spacer, p6, spacer, p7, spacer,
                  spacer, spacer, spacer, spacer, spacer,
                  align="hv", ncol=5, nrow=5, rel_widths=c(0.075, 1, 0.075, 1, 0.05), 
                  rel_heights=c(0.05, 1, 0.1, 1, 0.1), axis="rltb", scale=1) +
  # y-axis labels
  draw_label("Infected host density ± SE", x=0.03, y=0.75, angle=90, size=14) +
  draw_label("Total cercarial production ± SE", x=0.5, y=0.75, angle=90, size=14) +
  draw_label("Per capita cercarial production ± SE", x=0.03, y=0.25, angle=90, size=14) +
  draw_label("Per capita cercarial production ± SE", x=0.5, y=0.25, angle=90, size=14) +
  # x-axis labels
  draw_label("Weeks post snail introduction", x = 0.55, y = 0.51, vjust=0 ,size=14) +
  draw_label("Weeks post snail introduction", x = 0.27, y = 0.02, vjust=0 ,size=14) +
  draw_label("Snail biomass density, mg", x = 0.78, y = 0.02, vjust=0 ,size=14) +
  # panel labels
  draw_label("A", x = 0.12, y = 0.95) + 
  draw_label("B", x = 0.62, y = 0.95) +
  draw_label("C", x = 0.12, y = 0.48) +
  draw_label("D", x = 0.62, y = 0.48) +
  # fit statistics
  draw_label(expression(paste(R^2, " = 0.30")), x=0.24, y=0.93, size=12)+
  draw_label(expression(paste(R^2, " = 0.66")), x=0.7, y=0.93, size=12)+
  draw_label(expression(paste(R^2, " = 0.67")), x=0.24, y=0.47, size=12)+
  # legend
  draw_label(expression(underline("Founder size")), x=0.86, y=0.95, size=11) +
  draw_label("●", x=0.81, y=0.93, size=16, colour="black") + # unicode character, found using character map
  draw_label("▲", x=0.81, y=0.91, size=12, colour="black") + # unicode character
  draw_label("▪", x=0.81, y=0.89, size=24, colour="black") + # unicode character
  draw_line(x = c(0.825, 0.86), y=c(0.9275, 0.9275), linetype="dotted") + draw_label("Small", x = 0.87, y=0.93, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.9075, 0.9075), linetype="dashed") + draw_label("Medium", x = 0.87, y=0.91, size=11, hjust=0) +
  draw_line(x = c(0.825, 0.86), y=c(0.8875, 0.8875), linetype="solid") + draw_label("Large", x = 0.87, y=0.89, size=11, hjust=0) +
  draw_label(expression(underline("Nutrients")), x=0.86, y=0.86, size=11) +
  draw_label("High", x=0.86, y=0.84, size=11, colour="limegreen") +
  draw_label("Low", x=0.86, y=0.82, size=11) 

save_plot("Fig4.png", Fig3B, ncol=3, nrow=4, base_height=2, base_aspect_ratio = 1.1, dpi=300, units="in")


