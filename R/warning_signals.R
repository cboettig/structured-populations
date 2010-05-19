# warning_signals.R code
# Author: Carl Boettiger <cboettig@gmail.com>
# 
# in src/ directory, run R CMD SHLIB warning_signals.c gillespie.c odeintegrators.c
# 

#require(Kendall)
#dyn.load("../src/warning_signals.so")


# sample_time should be called window
warning_signals <- function(sample_time = 50, sample_freq = 1, max_time = 270, n_ensembles = 1, start_polluting = 200, pollute_timestep = 1, pollute_increment = 1){


	N <- as.integer(max_time/sample_freq)

	out <- .C(	"warning_signals", 
				double(N), 
				double(N), 
				double(N), 
				double(N), 
				double(N), 
				double(N), 
				double(N), 
				as.integer(sample_time), 
				as.double(sample_freq), 
				as.integer(max_time), 
				as.integer(n_ensembles), 
				as.double(start_polluting), 
				as.double(pollute_timestep), 
				as.double(pollute_increment) )

	times <- out[[1]]
	enviro <- out[[2]]
	means <- out[[3]]
	vars <- out[[4]]
	skews <- out[[5]]
	ar1 <- out[[6]]
	arN  <- out[[7]]

	x <- times != 0

	output <- data.frame(	times = times[x],
							enviro = enviro[x],
							means = means[x],
							vars = vars[x],
							skews = skews[x],
							ar1 = ar1[x],
							arN = arN[x] )
	class(output) = "warning_signals"
	output
}

plot.warning_signals <- function(output){ 
#	x <- output$times != 0
#	burnin <- which(output$times == 2*output$times[1])
#	l <- sum(x)
	par(mfrow=c(3,2) )

	test_m <- Kendall(output$times, output$means)
	test_v <- Kendall(output$times, output$vars)
	test_s <- Kendall(output$times, output$skews)
#	test_1 <- Kendall(output$times[burnin:l], output$ar1[burnin:l])
	test_1 <- Kendall(output$times, output$ar1)
	test_N <- Kendall(output$times, output$arN)
	plot(output$times, output$means, type='l', xlab="time", ylab="means", main=paste("2-tail P: ", test_m$sl))
	plot(output$times, output$vars, type='l', xlab="time", ylab="var", main=paste("2-tail P: ", test_v$sl))
	plot(output$times, output$skews, type='l', xlab="time", ylab="skew", main=paste("2-tail P: ", test_s$sl))
#	plot(output$times[burnin:l], output$ar1[burnin:l], type='l', xlab="time", ylab="ar1", main=paste("2-tail P: ", test_1$sl))
	plot(output$times, output$ar1, type='l', xlab="time", ylab="ar1", main=paste("2-tail P: ", test_1$sl))
	plot(output$times, output$arN, type='l', xlab="time", ylab="arN", main=paste("2-tail P: ", test_N$sl))
	plot(output$times, output$enviro, type='l', xlab="time", ylab="enviro")

	print(c("2-sided P-value:", test_1$sl))
#	print(var(output$means[x]) )
}


#output <- warning_signals(sample_time = 50, max_time=300, n_ensembles=1, start_polluting=200, pollute_timestep = .1, pollute_increment=.05)

#plot(output)
