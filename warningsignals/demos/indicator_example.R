# indicator_example.R

#source("load_CaCO3.R")
#source("bootstrap_indicators.R")

#m<-fit_models(X, "LTC")
require(warningsignals)
require(socialR)
load("5550699471.Rdat")

social_plot(bootstrap_indicators(X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", width=480, height=4*480)

