#lin_bifur_models.R
source("../R/likelihood_bifur_models.R")
source("../R/gaussian_process.R")
require(pmc)
require(odesolve)

## Simulate a dataset under slow linear change
pars <- c(Ro=1, m= -0.01, theta=1, sigma=1)
X <- simulateGauss(timedep_LTC, pars, N=500, T=10)
plot(X)

## fit both const and timedep models
start <- c(Ro=.5, m=0, theta=.1, sigma=.1)
timedep <- updateGauss(timedep_LTC, start, X, control=list(maxit=1000))
start <- c(Ro=.5, theta=.1, sigma=.1)
const <- updateGauss(const_LTC, start, X, control=list(maxit=1000))

#update(timedep, simulate(timedep))
#update(const, simulate(const))

out <- montecarlotest(const, timedep, cpu=16)
save(list=ls(), file="lin_bifur_models.Rdat")

png("lin_bifur_models.png")
plot(out)
dev.off()


gitcom <- system('git commit -a -m "autocommit warning signals"', intern=TRUE)[[1]]
system(paste('flickr_upload --tag="stochpop warningsignals" --description="', gitcom,  '" lin_bifur_models.png', sep=""))
system(paste('hpc-autotweets "@cboettig #stochpop warningsignals done ', gitcom, '"', sep=""))




## Now just need an example where R(t) is constant, requires redefining LTC and LSN

