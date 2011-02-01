#indicator_vs_likelihood.R
tags <- "warningsignals stochpop"
nboot <- 16
require(socialR)
require(warningsignals)
sfInit(parallel=TRUE, cpu=16)
nboot=160
sfLibrary(warningsignals)


pars <- c(Ro=5.0, m= -.01, theta=100, sigma=1)
sfExportAll()
## Simulate a dataset under slow linear change
warning <- simulateGauss(timedep_LSN, pars, N=500, T=100, Xo=100)
no_warning <- simulateGauss(const_LSN, pars, N=500, T=100, Xo=100)

# Likelihood Fits to each data-set and their relative model comparison  
timedep <- updateGauss(timedep_LSN, pars, warning, control=list(maxit=1000))
const <- updateGauss(const_LSN, pars, warning, control=list(maxit=1000))
llik_warning <- 2*(loglik(timedep)-loglik(const))

timedep_no <- updateGauss(timedep_LSN, pars, no_warning, control=list(maxit=1000))
const_no <- updateGauss(const_LSN, pars, no_warning, control=list(maxit=1000))
llik_nowarning <- 2*(loglik(timedep_no)-loglik(const_no))
save(list=ls(), file="LSN.Rdat")

out <- montecarlotest(const, timedep, cpu=16, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="LSN.Rdat")
social_plot(plot(out), file="LSN.png", tag="warningsignal stochpop LSN")

