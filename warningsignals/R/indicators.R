# indicators.R
# simple functions to calculate statistics often used as leading indicators of a regime shift

## Dependencies:
## requires psych package for skew & kurosis calculations
## used to require Kendall package, has been replaced with cor.test from the base package, 
## which resolves ties for more accurate p value, but doesn't matter since we focus on tau 
## and becayse cor.test is used by others i.e. Dakos 2008

require(psych)

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

window_ar.ols <- function(X, windowsize=length(X)/2, demean=FALSE){
	sapply(0:(length(X)-windowsize), 
		function(i){
			a<-ar.ols(X[(i+1):(i+windowsize)], demean=demean)
			a$ar[1]
		}
	)
}


warning_stats <- function(X, indicator){
	if(is(X,"ts")){
		w <- length(X)/2
		end <- length(X)
#		out <- Kendall(time(X)[w:end], indicator(X))
		out <- cor.test(time(X)[w:end], indicator(X), method="kendall")
	} else {
		w <- length(X[,1])/2
		end <- length(X[,1])
#		out <- Kendall(X[w:end,1], indicator(X[,2]))
		out <- cor.test(X[w:end,1], indicator(X[,2]), method="kendall")
	}
#	p <- as.numeric(out$sl)
#	c(as.numeric(out$tau), p)
	c(out$estimate, out$p.value)
}




