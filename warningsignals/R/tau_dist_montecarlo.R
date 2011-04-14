############## Define a bunch of useful plotting functions etc #############################
## Type I & II error rates for the distribution of the test statistic
err_rates <- function(null_dist, test_dist, p=.05){
	sig <- null_dist[2,] < p
	null_err <- sum(null_dist[1,sig] > 0)/length(null_dist[1,])
	sig <- test_dist[2,] < p
	test_err <- sum(test_dist[1,sig] < 0)/length(test_dist[1,])
	c(null_err=null_err, test_err=test_err)
}

plt_tau <- function(test_tau_dist, null_tau_dist, indicator, ylim=NULL, legend=FALSE, show_p=TRUE, ...){
	td <- density(test_tau_dist[1,])
	nd <- density(null_tau_dist[1,])
	if(is.null(ylim)) ylim <- c( min(nd$y, td$y), max(nd$y, td$y))
	plot(nd, type="n", main="", col=rgb(0,0,1,1), xlim=c(-1,1), ylim=ylim, ...)
	polygon(nd$x, nd$y, col=rgb(0,0,1,.3), border=rgb(0,0,1,.5))
	polygon(td$x, td$y, col=rgb(1,0,0,.3), border=rgb(1,0,0,.5))
#	polygon(nd$x, nd$y, col='blue', density=8, lwd=2, border=NA)
#	polygon(td$x, td$y, col='pink', density=8, lwd=2, border=NA, angle=-45)
#	lines(nd, lwd=2) #, col="lightblue")
#	lines(td, lwd=2) #, col="pink")

	#lines(nd, col="blue", lwd=3)
	if(legend) legend("topright", c("test", "null"), pch=15, col=c("red", "blue"))
	if(show_p){
		text(xshift(0), yshift(5), paste("fraction of test with p <0.05 is ", sum(test_tau_dist[2,] <.05)/length(null_tau_dist[2,])), pos=4 )
		text(xshift(0), yshift(15), paste("frac of null with p <0.05 is ", sum(null_tau_dist[2,] <.05)/length(null_tau_dist[2,])), pos=4 )
	}
}


tau_dist_montecarlo <- function(X, const, timedep, signal=c("Variance", "Autocorrelation", "Skew", "Kurtosis"), nboot=200, cpu=2, windowsize=round(length(X)/2), method=c("kendall", "pearson", "spearman"))
## Compute Monte Carlo bootstrap of tau under each model
{

	## display the simple likelihood comparison between models
	## be worried if this is very negative, or very large, could
	## indicate poor convergence
	print( llik_warning_fit <- 2*(loglik(timedep)-loglik(const)) )

	method = match.arg(method)
	signal = match.arg(signal)
	observed <- compute_tau(X, signal, windowsize, method=method) 

	## prepare parallel environment
	if(cpu>1 & !sfIsRunning()){ 	
		sfInit(parallel=TRUE, cpu=cpu)
		sfLibrary(warningsignals)
		sfExportAll()
	} else if(cpu<2 & !sfIsRunning()){
		sfInit()
	} 

## Look at the distribution of Taus when simulating from timedep
	test_tau_dist <- sfSapply(1:nboot, function(i){
		Z <- simulate(timedep)
		compute_tau(Z, signal, windowsize, method=method)
	})
## Distribution of Taus simulating from const model
	null_tau_dist <- sfSapply(1:nboot, function(i){
		Y <- simulate(const)
		compute_tau(Y, signal, windowsize, method=method)
	})

## should pass out a generic "observed"
	out <- list(test_tau_dist=test_tau_dist, null_tau_dist=null_tau_dist,
		signal=signal, X=X, llik_warning_fit=llik_warning_fit,
		const=const, timedep=timedep, observed = observed, windowsize=windowsize, method=method)
	class(out) <- "tau_dist_montecarlo"
	out
}

## should combine with plt_tau(?)
plot.tau_dist_montecarlo <- function(out, show_p=TRUE, show_error=FALSE, threshold=.95, ...){

	plt_tau(out$test_tau_dist, out$null_tau_dist, out$signal, show_p=show_p, ...)
#	abline(v=out$observed[1], lty=2, lwd=3, col="darkred")
	points(out$observed[1],yshift(1), cex=1.5, col="black", pch=25, fg="black", bg="black")
	if(show_error){	
		## Power/Error calculation
		test_dist <- unlist(out$test_tau_dist[1,])
		null_dist <-  unlist(out$null_tau_dist[1,])
		nboot <- length(test_dist)
		threshold_tail <- sort(null_dist)[ round(threshold*nboot) ]
		power <- sum(test_dist > threshold_tail)/nboot
		p <- 1-sum(null_dist < out$observed[1])/nboot
		text(xshift(0), yshift(95), paste("Type I:", round(p,3), "\n Type II:", round(1-power,3)), pos=4,
             cex=par()$cex.lab/.8)
	}
}























###################### Depricated plots ###

## Plot data for a single input set  
plot_kendalls <- function(X){

	if(!is.ts(X)){
		n <- length(X[,1])
		start <- X[1,1]
		end <- X[1,n]
		X <- ts(X[,2], start=start, end=end, freq=(end-start)/n)
	}

	par(mfrow=c(3,1))
#	plot(R(time(X), pars), col="red")

	plot(window_var(X), type="l", main="Variance", xlab="Time")
	show_stats(X, window_var)

	plot(window_autocorr(X), type="l", main="Autocorrelation", xlab="Time")
	show_stats(X, window_autocorr)
 
}





## Wrapper to plot_data, simulates a sample from each model
plot_sample <- function(const,timedep){
		## Creates some simulated data from these estimates and show example performance
		no_warning <- simulate(const)
		warning <- simulate(timedep)
		plt_data(warning, no_warning)
}

## Collective panel plots: Warning signals for the actual warning signal and the constant conditions model.  Should be modified to match plot.dakos format...
plt_data <- function(warning, no_warning){
## Description: plots two datasets, their variances and autocorrelations in a 3x2 frame
## Args: 
##			warning -- a dataset (ts or matrix)
##			no_warning -- another dataset (ts or matrix)
	par(mfrow=c(3,2))
	plot(warning, main="Stability loss (LTC model)")
	plot(no_warning, main="No stability loss (const eval)")
#	plot(R(time(X), pars), col="red")

	if(is(warning, "ts")){
		plot(window_var(warning), type="l", main="Variance", xlab="Time")
		show_stats(warning, window_var)

		plot(window_var(no_warning), type="l", main="Variance", xlab="Time")
		show_stats(no_warning, window_var)

		plot(window_autocorr(warning), type="l", main="Autocorrelation", xlab="Time")
		show_stats(warning, window_autocorr)

		plot(window_autocorr(no_warning), type="l", main="Autocorrelation", xlab="Time")
		show_stats(no_warning, window_autocorr)
	} else {
		plot(window_var(warning[,2]), type="l", main="Variance", xlab="Time")
		show_stats(warning, window_var)

		plot(window_var(no_warning[,2]), type="l", main="Variance", xlab="Time")
		show_stats(no_warning, window_var)

		plot(window_autocorr(warning[,2]), type="l", main="Autocorrelation", xlab="Time")
		show_stats(warning, window_autocorr)

		plot(window_autocorr(no_warning[,2]), type="l", main="Autocorrelation", xlab="Time")
		show_stats(no_warning, window_autocorr)
	} 
}


