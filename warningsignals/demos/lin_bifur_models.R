#lin_bifur_models.R
source("../R/likelihood_bifur_models.R")
source("../R/gaussian_process.R")
require(pmc)
require(odesolve)

## Simulate a dataset under slow linear change

n <- seq(10,200, by=10)
sfInit(parallel=TRUE, cpu=16)
data <- sfLapply(n, 

	function(i){
	pars <- c(Ro=50, m= -5*.9, theta=1, sigma=1)
	X <- simulateGauss(timedep_LTC, pars, N=N[i], T=10)

	## fit both const and timedep models
	start <- c(Ro=.5, m=0, theta=.1, sigma=.1)
	timedep <- updateGauss(timedep_LTC, start, X, control=list(maxit=1000))
	start <- c(Ro=.5, theta=.1, sigma=.1)
	const <- updateGauss(const_LTC, start, X, control=list(maxit=1000))

	out <- montecarlotest(const, timedep, cpu=16)
	save(list=ls(), file="lin_bifur_models.Rdat")
	png("timeseries.png")
	plot(X)
	dev.off()
	png("lin_bifur_models.png")
	plot(out)
	dev.off()

	id <- i 
	gitcom <- system('git log -n -1', intern=TRUE)[[1]]
	system(paste('flickr_upload --tag="stochpop warningsignals" --description="', gitcom,  '" lin_bifur_models.png timeseries.png', sep=""))
	system(paste('hpc-autotweets "#stochpop iteration id = ', id, gitcom, '"', sep=""))

	out
})

save(list=ls(), file="lin_bifur_models.Rdat")

png("n_power_curve.png")
plot(n, sapply(1:length(n), function(i) data[[i]]$power))
dev.off()
	gitcom <- system('git log -n -1', intern=TRUE)[[1]]
	system(paste('flickr_upload --tag="stochpop warningsignals" --description="', gitcom,  '" n_power_curve.png', sep=""))
	system(paste('hpc-autotweets "@cboettig #stochpop warningsignals done, id = ', id, gitcom, '"', sep=""))



## Now just need an example where R(t) is constant, requires redefining LTC and LSN

