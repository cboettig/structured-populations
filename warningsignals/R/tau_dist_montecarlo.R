############## Define a bunch of useful plotting functions etc #############################

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
## Type I & II error rates for the distribution of the test statistic
err_rates <- function(null_dist, test_dist, p=.05){
	sig <- null_dist[2,] < p
	null_err <- sum(null_dist[1,sig] > 0)/length(null_dist[1,])
	sig <- test_dist[2,] < p
	test_err <- sum(test_dist[1,sig] < 0)/length(test_dist[1,])
	c(null_err=null_err, test_err=test_err)
}

plt_tau <- function(test_tau_dist, null_tau_dist, type){
	plot(density(test_tau_dist[1,]), main=paste("Kendall's Tau in ", type), lwd=3, col="blue", xlim=c(-1,1))
	lines(density(null_tau_dist[1,]), col="red", lwd=3)
	legend("topright", c("test", "null"), lwd=3, col=c("blue", "red"))
	text(xshift(2), yshift(1.5), paste("frac test p <0.05 is ", sum(test_tau_dist[2,] <.05)/length(null_tau_dist[2,])), cex=1.5, font=2, col="blue")
	text(xshift(2), yshift(1), paste("frac null p <0.05 is ", sum(null_tau_dist[2,] <.05)/length(null_tau_dist[2,])), cex=1.5, font=2, col="red")
}


## Collective panel plots: Warning signals for the actual warning signal and the constant conditions model
plt_data <- function(warning, no_warning){
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


tau_dist_montecarlo <- function(X, const, timedep, signal=c("Variance", "Autocorrelation"), nboot=200, cpu=2){
	print( llik_warning_fit <- 2*(loglik(timedep)-loglik(const)) )

	sfInit(parallel=TRUE, cpu=cpu)
	sfLibrary(warningsignals)
	sfExportAll()
## Look at the distribution of Taus
	test_tau_dist <- sfSapply(1:nboot, function(i){
		X <- simulate(timedep)
		warning_stats(X, window_var)
	})

	null_tau_dist <- sfSapply(1:nboot, function(i){
		Y <- simulate(const)
		warning_stats(Y, window_var)
	})
	out <- list(test_tau_dist=test_tau_dist, null_tau_dist=null_tau_dist,
		signal=signal, X=X, llik_warning_fit=llik_warning_fit,
		const=const, timedep=timedep)
	class(out) <- "tau_dist_montecarlo"
	out
}

plot.tau_dist_montecarlo <- function(out, show_sample=FALSE){
	if(show_sample){ 
		## Createis some simulated data from these estimates and show example performance
		warning <- simulate(out$timedep)
		no_warning <- simulate(out$const)
		plt_data(warning, no_warning)
	} else {
		plt_tau(out$test_tau_dist, out$null_tau_dist, out$signal)
	}
}




