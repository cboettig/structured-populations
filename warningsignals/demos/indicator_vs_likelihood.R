#indicator_vs_likelihood.R
nboot <- 800
cpu <- 16
require(socialR)
require(warningsignals)
gitcommit()


########################### Begin actual analysis ######################## 
pars <- c(Ro=5.0, m= -.05, theta=100, sigma=1)
const_pars <- c(Ro=5.0, theta=100, sigma=1)


## Some initial data: Simulate some sample data under slow linear change 
X <- simulateGauss(timedep_LSN, pars, N=100, T=99, Xo=100)

## traditional stats on this

leading_indicators <- function(X){
## mgp is margin of title, axis label, axis line.  3,1,0 is default
	par(cex.lab=1.7, lwd=2, mgp=c(2,.4,0) )

#	par(mfrow=c(5,1))
## postitions of the plots 1, 2, 3, 4, 5 in a matrix layout
	mat <-	rbind(c(1),c(2), c(3), c(4), c(5) )
	layout(mat, height = c(1.4,1,1,1,1.45))

## mar is margins, in order bottom, left, top, right.  default is 5,4,4,2
	par( mar=c(0,6,4,2) ) ## top margin
	plot(X, type="o", xaxt="n", ylab="data")
	par( mar=c(0,6,0,2) ) ## no top or bottom margin
	plot_indicator(X, "Variance", xaxt="n")
	plot_indicator(X, "Autocor", xaxt="n")
	plot_indicator(X, "Skew", xaxt="n")
	par( mar=c(5,6,0,2) ) ## restore bottom margin
	plot_indicator(X, "Kurtosis")
}
leading_indicators(X)
#plot_indicator(X, "CV")



## grab some info for reporting
T <- max(time(X))
N <- length(X)
Xo <- X@.Data[1]
sampling <- c(T=T, N=N, Xo=Xo)

comment <- paste(c(names(pars), ":", pars, "\\", names(sampling), ":", sampling, "nboot:", nboot, " Transcritical"), collapse=" ")

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
out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
social_plot(plot(out), file="indicator_vs_likelihood_mc.png", tag="warningsignals stochpop MC LTC", comment=comment)



