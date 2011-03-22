require(socialR)
require(warningsignals)
gitcommit()


########################### Begin actual analysis ######################## 
pars <- c(Ro=5.0, m= -.04999, theta=500, sigma=5)
const_pars <- c(Ro=5.0, theta=500, sigma=5)
## Some initial data: Simulate some sample data under slow linear change 
X <- simulateGauss(timedep_LSN, pars, N=100, T=100, Xo=500)

all_indicators(X)

## grab some info for reporting
T <- max(time(X)); N <- length(X); Xo <- X@.Data[1]; sampling <- c(T=T, N=N, Xo=Xo)
comment <- paste(c(names(pars), ":", pars, "\\", names(sampling), ":", sampling, "nboot:", nboot, "LSN"), collapse=" ")
## traditional stats on this
social_plot(all_indicators(X), file="indicators.png", tags="warningsignals stochpop tau indicators", comment=comment)



