# indicators.R
# simple functions to calculate statistics often used as leading indicators of a regime shift

window_var <- function(X, windowsize=length(X)/2){
	sapply(0:(length(X)-windowsize), function(i){
		var(X[(i+1):(i+windowsize)]) 
	})
}


window_mean <- function(X, windowsize=length(X)/2){
	sapply(0:(length(X)-windowsize), function(i){
		mean(X[(i+1):(i+windowsize)]) 
	})
}

require(psych)
window_skew <- function(X, windowsize=length(X)/2){
	sapply(0:(length(X)-windowsize), function(i){
		skew(X[(i+1):(i+windowsize)]) 
	})
}

window_kurtosi <- function(X, windowsize=length(X)/2){
	sapply(0:(length(X)-windowsize), function(i){
		kurtosi(X[(i+1):(i+windowsize)]) 
	})
}


window_autocorr <- function(X, windowsize=length(X)/2){
	sapply(0:(length(X)-windowsize), 
		function(i){
			a<- acf(X[(i+1):(i+windowsize)], lag.max=1, plot=F) 
			a$acf[2]
		}
	)
}


warning_stats <- function(X, indicator){
	if(is(X,"ts")){
		w <- length(X)/2
		end <- length(X)
		out <- Kendall(time(X)[w:end], indicator(X))
	} else {
		w <- length(X[,1])/2
		end <- length(X[,1])
		out <- Kendall(X[w:end,1], indicator(X[,2]))
	}
#	p <- format.pval(out$sl)
	p <- as.numeric(out$sl)
	c(as.numeric(out$tau), p)
}




