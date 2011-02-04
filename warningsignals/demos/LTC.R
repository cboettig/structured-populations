#indicator_vs_likelihood.R
tags <- "warningsignals stochpop"
cpu <- 16
nboot <- 16
require(socialR)
require(warningsignals)
#sfInit(parallel=TRUE, cpu=16)
#sfLibrary(warningsignals)
pars <- c(Ro=5.0, m= -.04, theta=100, sigma=1)
#sfExportAll()
## Simulate a dataset under slow linear change
warning <- simulateGauss(timedep_LTC, pars, N=500, T=100, Xo=100)

# Likelihood Fits to each data-set and their relative model comparison  
timedep <- updateGauss(timedep_LTC, pars, warning, control=list(maxit=1000))
const <- updateGauss(const_LTC, pars, warning, control=list(maxit=1000))
llik_warning <- 2*(loglik(timedep)-loglik(const))


out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="LTC.Rdat")
social_plot(plot(out), file="LTC.png", tag="warningsignal stochpop LTC", comment=paste("pars = ", pars))

