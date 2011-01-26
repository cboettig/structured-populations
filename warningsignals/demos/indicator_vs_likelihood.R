#indicator_vs_likelihood.R
# edit stuff 
tags <- "warningsignals stochpop"
require(warningsignals)
require(socialR)
tweet_errors(tags=tags)

## Simulate a dataset under slow linear change

sfInit(parallel=TRUE, cpu=16)
sfLibrary(warningsignals)
sfLibrary(socialR)
sfExportAll()

pars <- c(Ro=50, m= -45, theta=1, sigma=1)
X <- simulateGauss(timedep_LTC, pars, N=500, T=1)
plot(X)



