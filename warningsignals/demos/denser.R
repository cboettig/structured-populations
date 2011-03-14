#denser.R
nboot <- 160
cpu <- 16
require(socialR)
gitcommit()
require(warningsignals)


########################### Begin actual analysis ######################## 
pars <- c(Ro=5.0, m= -.015, theta=100, sigma=1)
const_pars <- c(Ro=5.0, theta=100, sigma=1)
## Some initial data: Simulate some sample data under slow linear change 
X <- simulateGauss(timedep_LTC, pars, N=2500, T=100, Xo=100)

# initialize and fit models (could just use updateGauss instead of generic
#const <- updateGauss(const_LTC, const_pars, X, control=list(maxit=1000))  
const <- updateGauss(constOU, const_pars, X, control=list(maxit=1000))  
timedep <-updateGauss(timedep_LTC, pars, X, control=list(maxit=1000)) 

# Tau approach comparison
tau_var <- tau_dist_montecarlo(X, const, timedep, signal="Variance", nboot=nboot, cpu=cpu)
tau_acor <- tau_dist_montecarlo(X, const, timedep, signal="Autocorrelation", nboot=nboot, cpu=cpu)

save(list=ls(), file="denser.Rdat")
social_plot(plot(tau_var), file="taudist_var.png", tags="warningsignals stochpop tau var simulation", comment="N=2500, T=100, R=5-0.015*t", mention="cboettig")
social_plot(plot(tau_acor), file="taudist_acor.png", tags="warningsignals stochpop tau acor simulation", comment="N=2500, T=100, R=5-0.015*t")


## MONTECARLO Non-parametric bootstrap using the exact values.  
## As noted above, in real data we would use the MLEs, which has the point-estimate problem.  
out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="denser.Rdat")
social_plot(plot(out), file="denser_mc.png", tag="warningsignal stochpop LTC, simulation", comment="N=2500, T=100, R=5- 0/05*t")



