# indicator_example.R

fit_models <- function(X, method=c("LTC", "LSN")){
	const_pars <- c(Ro=1/max(time(X)), theta=mean(X), sigma=sd(X))
## Fit a linearized transcritical bifurcation model
	const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))
	pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
	timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))
	list(X=X, const=const, timedep=timedep, pars=pars, const_pars=const_pars)
}


bootstrap_indicators <- function(X, const, timedep, nboot=160, cpu=16){
# Tau approach comparison
	tau_var <- tau_dist_montecarlo(X, const, timedep, signal="Variance", nboot=nboot, cpu=cpu)
	tau_acor <- tau_dist_montecarlo(X, const, timedep, signal="Autocorrelation", nboot=nboot, cpu=cpu)
	tau_skew <- tau_dist_montecarlo(X, const, timedep, signal="Skew", nboot=nboot, cpu=cpu)
	tau_kurtosis <- tau_dist_montecarlo(X, const, timedep, signal="Kurtosis", nboot=nboot, cpu=cpu)

	par(mfrow=c(4,1))
	plot(tau_var)
	plot(tau_acor)
	plot(tau_skew)
	plot(tau_kurtosis)
}



