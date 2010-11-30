#lin_bifur_models.R
source("../R/likelihood_bifur_models.R")
source("../R/gaussian_process.R")

pars <- c(Ro=1, m=0, theta=1, sigma=1)
X <- simulate.gauss(setLTC, pars, N=500, T=10)
M <- update.gauss(setLTC, pars, X)


Y <- simulate.gauss(setLSN, pars)
