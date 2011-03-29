
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts
load("5555038554.Rdat")

source("../R/indicators.R")
indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
## FIGURE 1
social_plot(
	all_indicators(list(CaCO3=X, 
				   GlaciationI=data[[1]]$X_ts, 
				   GlaciationII=data[[2]]$X_ts ),	
				   indicators=indicators, method="kendall"),
	tags="warningsignals stochpop", 
	height=length(indicators)*240, width=3*480
)




