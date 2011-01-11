#lin_bifur_models.R
tags <- "warningsignals stochpop"


require(warningsignals)
require(socialR)

tweet_errors(tags=tags)

## Simulate a dataset under slow linear change

sfInit(parallel=TRUE, cpu=16)
sfLibrary(warningsignals)
sfLibrary(socialR)
sfExportAll()

	pars <- c(Ro=50, m= -.4, theta=1, sigma=1)
	X <- simulateGauss(timedep_LTC, pars, N=50, T=1)

	## fit both const and timedep models
	start <- c(Ro=.5, m=0, theta=.1, sigma=.1)
	timedep <- updateGauss(timedep_LTC, start, X, control=list(maxit=1000))
	start <- c(Ro=.5, theta=.1, sigma=.1)
	const <- updateGauss(const_LTC, start, X, control=list(maxit=1000))




	## Should use fixed parameters not those estimated from a simulation

#	const$pars <- c(Ro=.5, theta=1, sigma=1)
#	timedep$pars <- c(
	out <- montecarlotest(const, timedep, cpu=2, nboot=16)
#	save(list=ls(), file="lin_bifur_models.Rdat")

#	social_plot(plot(X), file="timeseries.png", tags=tags)
#	social_plot(plot(out), file="lin_bifur_models.png",tags=tags)


