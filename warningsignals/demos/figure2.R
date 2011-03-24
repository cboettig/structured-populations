## FIGURE 2
require(warningsignals)
require(socialR)
gitcommit()



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


## bootstrap each 
simLTC_taus <- bootstrap_tau(simLTC_m$X, simLTC_m$const, simLTC_m$timedep, indicators = c("Variance", "Autocorrelation"), nboot=2000, cpu=16)
CaCO3_taus <- bootstrap_tau(CaCO3_m$X, CaCO3_m$const, CaCO3_m$timedep, indicators = c("Variance", "Autocorrelation"), nboot=2000, cpu=16)
deut3_taus <- bootstrap_tau(deut3_m$X, deut3_m$const, deut3_m$timedep, indicators = c("Variance", "Autocorrelation"), nboot=2000, cpu=16)

save(list=ls(), file="figure2.Rdat")

social_plot(
	plot.bootstrap_tau(list(CaCO3=CaCO3_taus, Deut3=deut3_taus)),
	tags = "warningsingals stochpop", height=480*2, width=480*2
	)


