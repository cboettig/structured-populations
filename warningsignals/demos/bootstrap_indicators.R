# bootstrap_indicators.R

## A few wrapper functions to make it easy to bootstrap a secified set of indicator statistics

fit_models <- function(X, method=c("LTC", "LSN")){
	const_pars <- c(Ro=1/max(time(X)), theta=mean(X), sigma=sd(X))
## Fit a linearized transcritical bifurcation model
	const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))
	pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
	if(method=="LTC"){
		timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))
	} else if(method=="LSN"){
		timedep <- updateGauss(timedep_LSN, pars, X, control=list(maxit=1000))
	}
	list(X=X, const=const, timedep=timedep, pars=pars, const_pars=const_pars, method=method)
}


bootstrap_tau <- function(X, const, timedep, indicators = c("Variance", "Autocorrelation", "Skew", "Kurtosis"), nboot=160, cpu=16, windowsize=round(length(X)/2)){
# Tau approach comparison
	taus <- lapply(indicators, function(stat){ 	tau_dist_montecarlo(X, const, timedep, signal=stat, nboot=nboot, cpu=cpu) })
	class(taus) <- "boostrap_tau"
	
}


## This should really be able to plot the matrix of different models and different stats, with normalizing the height of the plots.
plot.bootstrap_tau <- function(taus){
	n <- length(taus)
	par(mfrow=c(n,1))
	for(i in 1:n) plot(taus[[i]])
}

