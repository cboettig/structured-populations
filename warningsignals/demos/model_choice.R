
choosemodel <- function(X, nboot=160, cpu=16){
const_pars <- c(Ro=1/max(time(X)), theta=mean(X), sigma=sd(X))
# Fit a linearized transcritical bifurcation model
	const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))
	pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
	LTC <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))
	LSN <- updateGauss(timedep_LSN, pars, X, control=list(maxit=1000))
	out <- montecarlotest(LTC, LSN, nboot=nboot, cpu=cpu)
	plot(out)
}

source("load_CaCO3.R")
social_plot(choosemodel(X), file="CaCO3_modelchoice.png", tag=paste(tags, "MC", "modelchoice")) 

source("load_deut.R")
social_plot(choosemodel(data[[3]]$X_ts), file="deut3_modelchoice.png", tag=paste("stochpop warningsignals deut", "MC", "modelchoice")) 



load("5550314677.Rdat")
social_plot(choosemodel(X), file="simdat_modelchoice.png", tag=paste(tags, "MC", "modelchoice")) 

