############## Define a bunch of useful plotting functions etc #############################
## Type I & II error rates for the distribution of the test statistic
err_rates <- function(null_dist, test_dist, p=.05){
	sig <- null_dist[2,] < p
	null_err <- sum(null_dist[1,sig] > 0)/length(null_dist[1,])
	sig <- test_dist[2,] < p
	test_err <- sum(test_dist[1,sig] < 0)/length(test_dist[1,])
	c(null_err=null_err, test_err=test_err)
}

plt_tau <- function(test_tau_dist, null_tau_dist, indicator){
	td <- density(test_tau_dist[1,])
	nd <- density(null_tau_dist[1,])
	ylim <- c( min(nd$y, td$y), max(nd$y, td$y))
	plot(nd, main=paste("Kendall's Tau in ", indicator), lwd=1, col=rgb(0,0,1,1), xlim=c(-1,1), ylim=ylim)
	polygon(nd$x, nd$y, col=rgb(0,0,1,.5), border=rgb(0,0,1,.5))
	polygon(td$x, td$y, col=rgb(1,0,0,.5), border=rgb(1,0,0,.5))

	#lines(nd, col="blue", lwd=3)
	legend("topright", c("test", "null"), pch=15, col=c("red", "blue"))
	text(xshift(5), yshift(5), paste("fraction of test with p <0.05 is ", sum(test_tau_dist[2,] <.05)/length(null_tau_dist[2,])), pos=4 )
	text(xshift(5), yshift(10), paste("frac of null with p <0.05 is ", sum(null_tau_dist[2,] <.05)/length(null_tau_dist[2,])), pos=4 )
}

tau_dist_montecarlo <- function(X, const, timedep, signal=c("Variance", "Autocorrelation", "Skew", "Kurtosis"), nboot=200, cpu=2){
	print( llik_warning_fit <- 2*(loglik(timedep)-loglik(const)) )

	observed_acor <- warning_stats(X, window_autocorr)
	observed_var <- warning_stats(X, window_var)
	observed_skew <- warning_stats(X, window_skew)
	observed_kurtosis <- warning_stats(X, window_kurtosi)

	if(cpu>1 & !sfIsRunning()){ 	
		sfInit(parallel=TRUE, cpu=cpu)
		sfLibrary(warningsignals)
		sfExportAll()
	} else if(cpu<2 & !sfIsRunning()){  sfInit()
	} else { }


## Look at the distribution of Taus
	test_tau_dist <- sfSapply(1:nboot, function(i){
		Z <- simulate(timedep)
		if(signal=="Variance"){  out <- warning_stats(Z, window_var)
		} else if (signal=="Autocorrelation") { out <- warning_stats(Z, window_autocorr) 
		} else if (signal =="Skew") { out <- warning_stats(Z, window_skew)
		} else if (signal =="Kurtosis") { out <- warning_stats(Z, window_kurtosi)
		} else { message("signal type not recognized")  }
		out
	})
	null_tau_dist <- sfSapply(1:nboot, function(i){
		Y <- simulate(const)
		if(signal=="Variance"){  out <- warning_stats(Y, window_var)
		} else if (signal=="Autocorrelation") { out <- warning_stats(Y, window_autocorr)
		} else if (signal =="Skew") { out <- warning_stats(Y, window_skew)
		} else if (signal =="Kurtosis") { out <- warning_stats(Y, window_kurtosi)

		} else { message("signal type not recognized")  }
		out
	})
	out <- list(test_tau_dist=test_tau_dist, null_tau_dist=null_tau_dist,
		signal=signal, X=X, llik_warning_fit=llik_warning_fit,
		const=const, timedep=timedep, observed_var = observed_var, observed_acor = observed_acor)
	class(out) <- "tau_dist_montecarlo"
	out
}

plot.tau_dist_montecarlo <- function(out){
		plt_tau(out$test_tau_dist, out$null_tau_dist, out$signal)
		if(out$signal == "Variance") observed_tau <- out$observed_var
		if(out$signal == "Autocorrelation") observed_tau <- out$observed_acor
		if(out$signal == "Skew") observed_tau <- out$observed_skew
		if(out$signal == "Kurtosis") observed_tau <- out$observed_kurtosis
		abline(v=observed_tau[1], lty=2, lwd=3, col="darkred")
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


