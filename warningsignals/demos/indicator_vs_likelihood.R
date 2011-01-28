#indicator_vs_likelihood.R
# edit stuff 
tags <- "warningsignals stochpop"
nboot <- 16
require(socialR)
require(warningsignals)
sfInit(parallel=TRUE, cpu=2)
sfLibrary(warningsignals)

pars <- c(Ro=5.0, m= -.049, theta=100, sigma=1)
sfExportAll()


## Simulate a dataset under slow linear change
warning <- simulateGauss(timedep_LTC, pars, N=500, T=100, Xo=100)
no_warning <- simulateGauss(const_LTC, pars, N=500, T=100, Xo=100)

# Likelihood Fits  
timedep <- updateGauss(timedep_LTC, pars=c(Ro=5, m=0, theta=100, sigma=1), warning, control=list(maxit=1000))
const <- updateGauss(const_LTC, pars, warning, control=list(maxit=1000))
llik_warning <- 2*(loglik(timedep)-loglik(const))

timedep_no <- updateGauss(timedep_LTC, pars, no_warning, control=list(maxit=1000))
const_no <- updateGauss(const_LTC, pars, no_warning, control=list(maxit=1000))
llik_nowarning <- 2*(loglik(timedep_no)-loglik(const_no))

out <- montecarlotest(const, timedep, cpu=16, nboot=nboot)
social_plot(plot(out), file="test.png", tag="warningsignal stochpop")

## a quick labling function
xshift <- function(xsteps){
	deltax <- (par()$xaxp[2]-par()$xaxp[1])/5
	par()$xaxp[1]+xsteps*deltax
}
yshift <- function(ysteps){
	deltay <- (par()$yaxp[2]-par()$yaxp[1])/5
	par()$yaxp[2]-ysteps*deltay
}
show_stats <- function(X, indicator){
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
social_plot(plts(), file="indicators.png", tags=tags)



## Look at the distribution of Taus
test_tau_dist <- sfSapply(1:nboot, function(i){
	X <- simulateGauss(timedep_LTC, pars, N=500, T=100, Xo=100)
	warning_stats(X, window_var)
})

null_tau_dist <- sfSapply(1:nboot, function(i){
	Y <- simulateGauss(const_LSN, pars, N=500, T=100, Xo=100)
	warning_stats(Y, window_var)
})

plt <- function(){
	plot(density(test_tau_dist[1,]), main="Kendall's Tau with (test) and without (null) destablizing", lwd=3, col="blue", xlim=c(-1,1))
	lines(density(null_tau_dist[1,]), col="red", lwd=3)
	legend("topright", c("test", "null"), lwd=3, col=c("blue", "red"))
	text(xshift(1), yshift(1.5), paste("frac test p <0.05 is ", sum(test_tau_dist[2,] <.05)/length(null_tau_dist[2,])), cex=1.5, font=2, col="blue")
	text(xshift(1), yshift(1), paste("frac null p <0.05 is ", sum(null_tau_dist[2,] <.05)/length(null_tau_dist[2,])), cex=1.5, font=2, col="red")
}
plt()
social_plot(plt(), file="taudist.png", tags=tags)




