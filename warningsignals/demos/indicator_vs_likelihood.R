#indicator_vs_likelihood.R
# edit stuff 
tags <- "warningsignals stochpop"
require(warningsignals)

## Simulate a dataset under slow linear change
pars <- c(Ro=50, m= -45, theta=1, sigma=1)
warning <- simulateGauss(timedep_LTC, pars, N=500, T=1)
no_warning <- simulateGauss(timedep_LTC, pars, N=500, T=1)


## a quick labling function
show_stats <- function(X, indicator){
	xshift <- function(xsteps){
		deltax <- (par()$xaxp[2]-par()$xaxp[1])/5
		par()$xaxp[1]+xsteps*deltax
	}
	yshift <- function(ysteps){
		deltay <- (par()$yaxp[2]-par()$yaxp[1])/5
		par()$yaxp[2]-ysteps*deltay
	}
	w <- warning_stats(X, indicator)
	text(xshift(1), yshift(1), paste("tau = ", round(w[1],2)), cex=1.5, col="red", font=2)
	text(xshift(1), yshift(2), paste("p = ", format.pval(w[2])), cex=1.5, col="blue", font=2)
}

## Collective panel plots
plts <- function(){
	par(mfrow=c(3,2))
	plot(warning, main="Stability loss (LTC model)")
	plot(no_warning, main="No stability loss (const eval)")
#	plot(R(time(X), pars), col="red")
	plot(window_var(warning), type="l", main="Variance", xlab="Time")
	show_stats(warning, window_var)

	plot(window_var(no_warning), type="l", main="Variance", xlab="Time")
	show_stats(no_warning, window_var)

	plot(window_autocorr(warning), type="l", main="Autocorrelation", xlab="Time")
	show_stats(warning, window_autocorr)

	plot(window_autocorr(no_warning), type="l", main="Autocorrelation", xlab="Time")
	show_stats(no_warning, window_autocorr)

}
plts()


test_tau_dist <- sapply(1:nboot, function(i){
	X <- simulateGauss(timedep_LSN, pars, N=500, T=1, Xo=6)
	warning_stats(X, window_var)
})

null_tau_dist <- sapply(1:nboot, function(i){
	Y <- simulateGauss(const_LSN, pars, N=500, T=1, Xo=6)
	warning_stats(Y, window_var)
})



plot(density(test_tau_dist[1,]), main="Kendall's Tau wtih and without warning")




