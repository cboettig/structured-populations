require(warningsignals)
require(socialR)
load("5550863240.Rdat")
m<-fit_models(X, "LTC")
taus <- bootstrap_tau(m$X, m$const, m$timedep)
social_plot(plot(taus), file="bootstrap_tau.png", height=480*4, tags=paste(tags, "sim", "LTC"))


