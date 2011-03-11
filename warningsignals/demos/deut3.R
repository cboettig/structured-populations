require(socialR)
require(warningsignals)
tags<-"warningsignals stochpop climatedata deut"
cpu <- 16
nboot <- 16

load("deut_data.Rdat")

X <- data[[3]]$X_ts

# consider smarter estimate for Ro
const_pars <- c(Ro=1/max(time(X)), theta=mean(X), sigma=sd(X))

## Fit a linearized transcritical bifurcation model
const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))
pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=as.numeric(const$pars["theta"]), sigma=as.numeric(const$pars["sigma"]))

out <- heteroOU(X,pars)
fit <- updateGauss("heteroOU", pars, X, control=list(maxit=1000))


timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))

## Fit the linearized saddle-node bifurcation model
#const <- updateGauss(const_LSN, const_pars, X, control=list(maxit=1000))
#smart estimates for pars
#pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
#timedep <- updateGauss(timedep_LSN, pars, X, control=list(maxit=1000))

print(llik_warning <- 2*(loglik(timedep)-loglik(const)))
sfInit(parallel=TRUE, cpu=cpu)
sfLibrary(warningsignals)
sfExportAll()

# Tau approach comparison
tau_var <- tau_dist_montecarlo(X, const, timedep, signal="Variance", nboot=nboot, cpu=cpu)
tau_acor <- tau_dist_montecarlo(X, const, timedep, signal="Autocorrelation", nboot=nboot, cpu=cpu)

save(list=ls(), file="deut3.Rdat")
social_plot(plot(tau_var), file="taudist_var.png", tags=paste(tags, "tau var"), mention="cboettig")
social_plot(plot(tau_acor), file="taudist_acor.png", tags=paste(tags, "tau acor"))

# plot example data
social_plot(plot(tau_var, show_sample=TRUE), tags="warningsignals stochpop tau deut3")

## MonteCarlo Cox's delta approach
out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="deut3.Rdat")
social_plot(plot(out), file="deut3.png", tag="warningsignals stochpop climatedata deut3", mention="cboettig")


