tags="warningsignals stochpop CaCO3"
cpu <- 2
nboot <- 160
require(warningsignals)
require(socialR)

# CaCO3 record
caco3 <- read.table("../data/caco3.txt")
# labels, time is in millions of years Before Present, so make negative
caco3 <- data.frame("MYYrs"=-caco3$V1, "CaCO3"=caco3$V2)

## window the data
g_ca <- caco3$MYYrs >= -39 & caco3$MYYrs <= -32  # Data with collapse (for plot)
p_ca <- caco3$MYYrs >= -39 & caco3$MYYrs < -34  # Data used in warning signal
X <- data.frame("time"=caco3$MYYrs[p_ca], "data"=caco3$CaCO3[p_ca])
# Rather annoying to have time backwards and negative, lets reverse this.
X <- data.frame("time"=rev(X[,1] - min(X[,1])), "data"=rev(X[,2]))


dat <- dakos_data_processing(X)
social_plot(plot.dakos(dat), file="plot_dakos.png", tag=tags)


X <- dat$X_ts
const_pars <- c(Ro=5.0, theta=mean(X), sigma=sd(X)*5*2)
## Fit a linearized transcritical bifurcation model
const <- updateGauss(const_LTC, const_pars, X, control=list(maxit=1000))
pars <- c(Ro=as.numeric(const$pars["Ro"]), m= as.numeric(const$pars["Ro"]*max(time(X))), theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))

## Fit the linearized saddle-node bifurcation model
#const <- updateGauss(const_LSN, const_pars, X, control=list(maxit=1000))
#pars <- c(Ro=as.numeric(const$pars["Ro"]), m= as.numeric(const$pars["Ro"]*max(time(X))), theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
#timedep <- updateGauss(timedep_LSN, pars, X, control=list(maxit=1000))

print(llik_warning <- 2*(loglik(timedep)-loglik(const)))
sfInit(parallel=TRUE, cpu=cpu)
sfLibrary(warningsignals)
sfExportAll()

# Tau approach comparison
tau_var <- tau_dist_montecarlo(X, const, timedep, signal="Variance", nboot=nboot, cpu=cpu)
tau_acor <- tau_dist_montecarlo(X, const, timedep, signal="Autocorrelation", nboot=nboot, cpu=cpu)

save(list=ls(), file="CaCO3.Rdat")
social_plot(plot(tau_var), file="taudist_var.png", tags=paste(tags, "tau var"), mention="cboettig")
social_plot(plot(tau_acor), file="taudist_acor.png", tags=paste(tags, "tau acor"), comment="Using interpolated and detrended data")

# plot example data
social_plot(plot(tau_var, show_sample=TRUE), tags="warningsignals stochpop tau CaCO3")

## MonteCarlo Cox's delta approach
out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="CaCO3.Rdat")
social_plot(plot(out), file="LTC_CaCO3.png", tag="warningsignals stochpop LTC climatedata CaCO3", mention="cboettig", comment="Using interpolated and detrended data")





## Falsely pretend it has even time intervals
#myyrs <- caco3$MYYrs[dat]
#end=myyrs[1]; start=myyrs[length(myyrs)]
#X <- ts(caco3[dat,2],start=start,end=end, freq=length(dat)/(end-start))
#plot(caco3[dat,1], caco3[dat,2])


