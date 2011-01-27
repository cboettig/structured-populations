require(socialR)
require(warningsignals)
sfInit(parallel=TRUE, cpu=16)
sfLibrary(warningsignals)
sfExportAll()

pars <- c(Ro=1, m= -0.45, theta=1, sigma=1)
X <- simulateGauss(timedep_LTC, pars, N=500, T=10)
 
## fit both const and timedep models
start <- c(Ro=.5, m=0, theta=.1, sigma=.1)
timedep <- updateGauss(timedep_LTC, start, X, control=list(maxit=1000))
start <- c(Ro=.5, theta=.1, sigma=.1)
const <- updateGauss(const_LTC, start, X, control=list(maxit=1000))
 
 
out <- montecarlotest(const, timedep, cpu=16, nboot=160)
save(file="test.Rdat", list=ls())
social_plot(plot(out), file="test.png", tag="warningsignal stochpop")
