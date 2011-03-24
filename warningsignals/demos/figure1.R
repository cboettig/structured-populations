require(warningsignals)
require(socialR)
gitcommit()

source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts
load("5555038554.Rdat")

## FIGURE 1
social_plot(
	all_indicators(list(deteriorating=deteriorating, constant=constant, deuterium=deut3), indicators = c("Variance", "Autocorrelation") ),
	tags="warningsignals stochpop"
	)


