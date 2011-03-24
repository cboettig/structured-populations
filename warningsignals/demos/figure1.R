source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts

## FIGURE 1
social_plot(
	all_indicators( data.frame(CaCO3, deut3) ),
	tags="warningsignals stochpop"
	)


