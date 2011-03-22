# indicator_example.R
# Run with a selected set of data from indicator_vs_likelihood.R simulation
load("5550815238.Rdat")

nboot <- 2000
cpu <- 16
require(socialR)
require(warningsignals)
gitcommit()

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



