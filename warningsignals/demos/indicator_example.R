# indicator_example.R

source("load_CaCO3.R")
source("bootstrap_indicators.R")

m<-fit_models(X, "LTC")
social_plot(bootstrap_indicators(X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png")

