# indicator_example.R

source("bootstrap_indicators.R")
source("load_deut.R")
m<-fit_models(data[[3]]$X_ts, "LTC")
social_plot(bootstrap_indicators(m$X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", height=480*4, tags=paste(tags, "deut3", "LTC"))

source("load_CaCO3.R")
m<-fit_models(X, "LTC")
social_plot(bootstrap_indicators(m$X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", height=480*4, tags=paste(tags, "CaCO3", "LTC"))

load("5550863240.Rdat")
m<-fit_models(X, "LTC")
social_plot(bootstrap_indicators(m$X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", height=480*4, tags=paste(tags, "sim", "LTC"))


## Try with LSN on data
source("load_deut.R")
m<-fit_models(data[[3]]$X_ts, "LSN")
social_plot(bootstrap_indicators(m$X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", height=480*4, tags=paste(tags, "deut3", "LSN"))


source("load_CaCO3.R")
m<-fit_models(X, "LSN")
social_plot(bootstrap_indicators(m$X, m$const, m$timedep, nboot=1000), file="bootstrap_tau.png", height=480*4, tags=paste(tags, "CaCO3", "LTC"))


