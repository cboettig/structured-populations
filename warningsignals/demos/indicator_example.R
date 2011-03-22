# indicator_example.R

source("load_CaCO3.R")
source("bootstrap_indicators.R")
social_plot(bootstrap_indicators(X), file="bootstrap_tau.png")

