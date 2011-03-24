## FIGURE 3

## load libraries
require(warningsignals)
require(socialR)
gitcommit()

## load datafiles
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts
load("5555038554.Rdat") # an LSN sim, stable and collapsing
#load("5550863240.Rdat") # an LTC sim, just collapsing
#simLTC <- X


## fit models
deterior_m<-fit_models(deteriorating, "LSN")
constant_m<-fit_models(constant, "LSN")
CaCO3_m <- fit_models(CaCO3, "LSN")
deut3_m <- fit_models(deut3, "LSN")



## montecarlo likelihood
cpu <- 8
nboot <- 2000

deterior_mc <- 
		montecarlotest(deterior_m$const, deterior_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)

constant_mc <- 
		montecarlotest(constant_m$const, constant_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)

CaCO3_mc <- 
		montecarlotest(CaCO3_m$const, CaCO3_m$timedep,
		cpu=cpu, nboot=nboot, GetParNames=FALSE)

deut3_mc <- 
		montecarlotest(deut3_m$const, deut3_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)


save(list=ls(), file="figure3.Rdat")

social_plot(plot(constant_mc), file="CaCO3.png", tag="warningsignals stochpop MC sim constant")
social_plot(plot(deterior_mc), file="CaCO3.png", tag="warningsignals stochpop MC sim deterior")
#social_plot(plot(simLTC_mc), file="simLTC.png", tag="warningsignals stochpop MC simLTC")
social_plot(plot(CaCO3_mc), file="CaCO3.png", tag="warningsignals stochpop MC CaCO3")
social_plot(plot(deut3_mc), file="Deut3.png", tag="warningsignals stochpop MC deut3")





