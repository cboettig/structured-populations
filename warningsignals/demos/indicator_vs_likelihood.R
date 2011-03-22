#indicator_vs_likelihood.R
nboot <- 800
cpu <- 16
require(socialR)
require(warningsignals)
gitcommit()


########################### Begin actual analysis ######################## 
pars <- c(Ro=5.0, m= -.004999, theta=500, sigma=5)
const_pars <- c(Ro=5.0, theta=500, sigma=5)
## Some initial data: Simulate some sample data under slow linear change 
X <- simulateGauss(timedep_LSN, pars, N=100, T=100, Xo=500)

## grab some info for reporting
T <- max(time(X)); N <- length(X); Xo <- X@.Data[1]; sampling <- c(T=T, N=N, Xo=Xo)
comment <- paste(c(names(pars), ":", pars, "\\", names(sampling), ":", sampling, "nboot:", nboot, "LSN"), collapse=" ")
all_indicators(X)

## traditional stats on this
social_plot(all_indicators(X), file="indicators.png", tags="warningsignals stochpop tau indicators", comment=comment)


# initialize and fit models (could just use updateGauss instead of generic
const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))  
timedep <-updateGauss(timedep_LTC, pars, X, control=list(maxit=1000)) 

# Tau approach comparison
tau_var <- tau_dist_montecarlo(X, const, timedep, signal="Variance", nboot=nboot, cpu=cpu)
tau_acor <- tau_dist_montecarlo(X, const, timedep, signal="Autocorrelation", nboot=nboot, cpu=cpu)
tau_skew <- tau_dist_montecarlo(X, const, timedep, signal="Skew", nboot=nboot, cpu=cpu)
tau_kurtosis <- tau_dist_montecarlo(X, const, timedep, signal="Kurtosis", nboot=nboot, cpu=cpu)

social_plot(plot(tau_var), file="taudist_var.png", tags="warningsignals stochpop tau var", comment=comment)
social_plot(plot(tau_acor), file="taudist_acor.png", tags="warningsignals stochpop tau acor", comment=comment)
social_plot(plot(tau_skew), file="taudist_skew.png", tags="warningsignals stochpop tau skew", comment=comment)
social_plot(plot(tau_kurtosis), file="taudist_kurtosis.png", tags="warningsignals stochpop tau kurtosis", comment=comment)


## MONTECARLO Non-parametric bootstrap using the exact values.  
## As noted above, in real data we would use the MLEs, which has the point-estimate problem.  
#out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
#social_plot(plot(out), file="indicator_vs_likelihood_mc.png", tag="warningsignals stochpop MC LSN", comment=comment)



