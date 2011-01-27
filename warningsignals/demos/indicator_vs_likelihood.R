#indicator_vs_likelihood.R
# edit stuff 
tags <- "warningsignals stochpop"
require(warningsignals)

## Simulate a dataset under slow linear change
pars <- c(Ro=50, m= -45, theta=1, sigma=1)
warning <- simulateGauss(timedep_LSN, pars, N=500, T=1, Xo=6)
no_warning <- simulateGauss(timedep_LSN, pars, N=500, T=1, Xo=6)

w <- warning_stats(warning, window_var)

plts <- function(){
	par(mfrow=c(2,2))
	plot(warning)
	plot(window_var(warning), main="running variance")
	text(50,.9*max(window_var(warning)), paste("tau = ", round(w[1],2)))
text(50,.8*max(window_var(warning)), paste("p = ", format.pval(w[2])))
	plot(no_warning)
	plot(window_var(no_warning), main="running variance")
}



test_tau_dist <- sapply(1:nboot, function(i){
	X <- simulateGauss(timedep_LSN, pars, N=500, T=1, Xo=6)
	warning_stats(X, window_var)
})

null_tau_dist <- sapply(1:nboot, function(i){
	Y <- simulateGauss(const_LSN, pars, N=500, T=1, Xo=6)
	warning_stats(Y, window_var)
})



plot(density(test_tau_dist[1,]), main="Kendall's Tau wtih and without warning")




