require(warningsignals)
require(socialR)
#load("5550863240.Rdat")
#source("bootstrap_indicators.R")
#m<-fit_models(X, "LTC")
#taus <- bootstrap_tau(m$X, m$const, m$timedep)
#save(list=ls(), file="test_taus.Rdat")
load("5553958875.Rdat")
#social_plot(plot(taus[[1]]), file="boot_tau.png", tags="test")
social_plot(plot(taus), file="bootstrap_tau.png", tags="test")


