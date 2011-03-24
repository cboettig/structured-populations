## FIGURE 2
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

indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
nboot <- 2000
cpu <- 16


## fit models
deterior_m<-fit_models(deteriorating, "LSN")
constant_m<-fit_models(constant, "LSN")
CaCO3_m <- fit_models(CaCO3, "LSN")
deut3_m <- fit_models(deut3, "LSN")


## bootstrap each 
constant_taus <- bootstrap_tau(constant_m$X, constant_m$const,
							constant_m$timedep,
							indicators = indicators,
							nboot=nboot, cpu=cpu)

deterior_taus <- bootstrap_tau(deterior_m$X,
							deterior_m$const, deterior_m$timedep,
							indicators = indicators,
							nboot=nboot, cpu=cpu)

CaCO3_taus <- bootstrap_tau(CaCO3_m$X, CaCO3_m$const, CaCO3_m$timedep,
							indicators = indicators,
							nboot=nboot, cpu=cpu)

deut3_taus <- bootstrap_tau(deut3_m$X, deut3_m$const, deut3_m$timedep,
							indicators = indicators,
							nboot=nboot, cpu=cpu)

social_plot(
	plot.bootstrap_tau(list(deteriorating=deterior_taus,
			constant=constant_taus, Deut3=deut3_taus), cex.axis=2),
			tags = "warningsingals stochpop", 
			height=480*2, width=480*3
	)


