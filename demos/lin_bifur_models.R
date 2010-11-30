#lin_bifur_models.R
source("../R/likelihood_bifur_models.R")
source("../R/gaussian_process.R")

gaussmodel <- setLTC
#gaussmodel <- setLSN

pars <- c(Ro=1, m=0, theta=1, sigma=1)
X <- simulate.gauss(gaussmodel, pars, N=500, T=10)
start <- c(Ro=.5, m=0, theta=.1, sigma=.1)
M <- update.gauss(gaussmodel, start, X, control=list(maxit=1000))
M$par

