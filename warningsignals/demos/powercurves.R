#lin_bifur_models.R
require(warningsignals)
require(socialR)
tags <- "warningsignals stochpop powercurves.R"
tweet_errors(tags=tags)

## Simulate a dataset under slow linear change

## n is the number of sample points
#n <- seq(10,200, by=10)

## M is a parameter for the rate of change in stability loss
M <- seq(-4.9, 0, by=.5)

sfInit(parallel=TRUE, cpu=16)
sfLibrary(warningsignals)
sfExportAll()

data <- lapply(M,

	function(i){
	pars <- c(Ro=50, m= i, theta=1, sigma=1)
	X <- simulateGauss(timedep_LTC, pars, N=100, T=10)

	## fit both const and timedep models
	start <- c(Ro=.5, m=0, theta=.1, sigma=.1)
	timedep <- updateGauss(timedep_LTC, start, X, control=list(maxit=1000))
	start <- c(Ro=.5, theta=.1, sigma=.1)
	const <- updateGauss(const_LTC, start, X, control=list(maxit=1000))

	out <- montecarlotest(const, timedep, cpu=16)
	save(list=ls(), file="powercurves.Rdat")
	social_plot(plot(out), file="powercurves.png", tags=tags, comment=paste("M = ", i))

	out
})

save(list=ls(), file="powercurves.Rdat")

social_plot(plot(M, sapply(1:length(M), function(i) data[[i]]$power)), file="powercurves.png")



