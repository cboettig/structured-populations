require(warningsignals)
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts
load("5550863240.Rdat") # consider replacing with a LSN sim
sim_LTC <- X

## FIGURE 1
social_plot(
	all_indicators( list(sim_LTC, CaCO3, deut3), indicators = c("Variance", "Autocorrelation") ),
	tags="warningsignals stochpop"
	)


