require(socialR)
require(warningsignals)
gitcommit()


pars <- c(Ro=5.0, m= -.04999, theta=500, sigma=5)
const_pars <- c(Ro=5.0, theta=500, sigma=5)

## Some initial data: Simulate some sample data under slow linear change 
deteriorating <- simulateGauss(timedep_LSN, pars, N=100, T=100, Xo=500)

## estimate pars from the simulated data, for as close a sim as possible
const_pars <- updateGauss(constOU, const_pars, deteriorating, control=list(maxit=1000))$pars
constant <- simulateGauss(constOU, const_pars, N=100, T=100, Xo=500)

# put both datasets into a list and plot side-by-side
dat <- list(deteriorating=deteriorating, constant=constant)
all_indicators(dat)

social_plot(all_indicators(dat), tags="stochpop warningsignals sim")




#T <- max(time(X)); N <- length(X); Xo <- X@.Data[1]; sampling <- c(T=T, N=N, Xo=Xo)
#comment <- paste(c(names(pars), ":", pars, "\\", names(sampling), ":", sampling, "nboot:", nboot, "LSN"), collapse=" ")
#social_plot(all_indicators(X), file="indicators.png", tags="warningsignals stochpop tau indicators", comment=comment)



