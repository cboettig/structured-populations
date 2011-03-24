## load datafiles
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts
load("5550863240.Rdat") # consider replacing with a LSN sim
simLTC <- X

## fit models
simLTC_m<-fit_models(simLTC, "LTC")
CaCO3_m <- fit_models(CaCO3, "LSN")
deut3_m <- fit_models(deut3, "LSN")


## montecarlo likelihood
cpu <- 16
nboot <- 160
simLTC_mc <- 
montecarlotest(simLTC_m$const, simLTC_m$timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)

CaCO3_mc <- 
montecarlotest(CaCO3_m$const, CaCO3_m$timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)

deut3_mc <- 
montecarlotest(deut3_m$const, deut3_m$timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)


save(list=ls(), file="figure3.Rdat")

social_plot(plot(simLTC_mc), file="simLTC.png", tag="warningsignals stochpop MC simLTC")
social_plot(plot(CaCO3_mc), file="CaCO3.png", tag="warningsignals stochpop MC CaCO3")
social_plot(plot(deut3_mc), file="Deut3.png", tag="warningsignals stochpop MC deut3")





