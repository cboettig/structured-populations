# indicator_example.R

source("bootstrap_indicators.R")
source("load_deut.R")
m<-fit_models(data[[3]]$X_ts, "LTC")
social_plot(bootstrap_indicators(m$X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", height=480*4)


