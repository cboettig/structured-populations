source("load_CaCO3.R")

cpu <- 2
nboot <- 160

const_pars <- c(Ro=5.0, theta=mean(X), sigma=sd(X)*5*2)
## Fit a linearized transcritical bifurcation model
const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))
pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))

## Fit the linearized saddle-node bifurcation model
#const <- updateGauss(const_LSN, const_pars, X, control=list(maxit=1000))
#pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
#timedep <- updateGauss(timedep_LSN, pars, X, control=list(maxit=1000))

print(llik_warning <- 2*(loglik(timedep)-loglik(const)))
sfInit(parallel=TRUE, cpu=cpu)
sfLibrary(warningsignals)
sfExportAll()

# Tau approach comparison
tau_var <- tau_dist_montecarlo(X, const, timedep, signal="Variance", nboot=nboot, cpu=cpu)
tau_acor <- tau_dist_montecarlo(X, const, timedep, signal="Autocorrelation", nboot=nboot, cpu=cpu)

social_plot(plot(tau_var), file="taudist_var.png", tags=paste(tags, "tau var"), mention="cboettig")
social_plot(plot(tau_acor), file="taudist_acor.png", tags=paste(tags, "tau acor"), comment="Using interpolated and detrended data")

# plot example data
social_plot(plot(tau_var, show_sample=TRUE), tags="warningsignals stochpop tau CaCO3")

## MonteCarlo Cox's delta approach
out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
social_plot(plot(out), file="LTC_CaCO3.png", tag="warningsignals stochpop LTC climatedata CaCO3", mention="cboettig", comment="Using interpolated and detrended data")





## Falsely pretend it has even time intervals
#myyrs <- caco3$MYYrs[dat]
#end=myyrs[1]; start=myyrs[length(myyrs)]
#X <- ts(caco3[dat,2],start=start,end=end, freq=length(dat)/(end-start))
#plot(caco3[dat,1], caco3[dat,2])


