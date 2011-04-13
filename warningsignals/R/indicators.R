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

window_cv <- function(X, windowsize=length(X)/2){
	sapply(0:(length(X)-windowsize), function(i){
			var(X[(i+1):(i+windowsize)]) / mean(X[(i+1):(i+windowsize)]) 
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


compute_indicator <- function(X, indicator=c("Autocorrelation", "Variance", "Skew", "Kurtosis", "CV"), windowsize=round(length(X)/2))
## Description Wrapper function to choose warning signal
## Assumes X is a ts object
{
	indicator = match.arg(indicator)
	if(indicator == "Autocorrelation"){
		out <- window_autocorr(X, windowsize)
	} else if(indicator == "Variance"){ 
		out <- window_var(X, windowsize)
	} else if(indicator == "Skew"){
		out <- window_skew(X, windowsize)
	} else if(indicator == "Kurtosis"){
		out <- window_kurtosi(X, windowsize)
	} else if(indicator == "CV"){
		out <- window_cv(X, windowsize)
	} else { warning(paste("Indicator", indicator, "not recognized")) 
	}
	out
}



## Compute and plot the given indicator
plot_indicator <- function(X, indicator=c("Autocorrelation", "Variance", "Skew", "Kurtosis", "CV"), windowsize=length(X)/2, xpos=0, ypos=90, method=c("kendall", "pearson", "spearman"), ...)
## Description
## Args:
##		X -- data, either ts object or matrix w/ time in col 1 and data in col 2
##		indicator -- name of the early warning indicator to be computed
##		windowsize -- size of the sliding window over which statistic is calculated
##		xpos, ypos -- % position from bottom left for text to appear
##		... -- extra arguments passed to plot
## Returns:
##		A plot of indicator statistic over time interval, with kendall's tau and p-val
{	
	method <- match.arg(method)
	if(!is.ts(X)){
		n <- length(X[,1])
		start <- X[1,1]
		end <- X[1,n]
		X <- ts(X[,2], start=start, end=end, freq=(end-start)/n)
	}
	Y <- compute_indicator(X, indicator, windowsize)
	time_window <- time(X)[windowsize:length(X)]
	plot(time_window, Y, xlim=c(start(X)[1], end(X)[1]), type="l", xlab="time", ylab=indicator, ...)
	abline(v=time_window[1], lty="dashed")

	out <- cor.test(time(X)[windowsize:length(X)], Y, method=method)

	w <- c(out$estimate, out$p.value)
    if (method=="kendall"){ 
## newline character doesn't work
#            text(xshift(xpos), yshift(ypos),
#                 substitute(paste(tau == val, "\n (p== ", pval, ")"), list(val=round(w[1],2),pval=format.pval(w[2]))), 
#            pos=4, cex=par()$cex.lab)
            text(xshift(0), yshift(ypos),
                 substitute(paste(tau == val), 
                            list(val=round(w[1],2) )),
                 pos=4, cex=par()$cex.axis)
            text(xshift(0), yshift(ypos-20),
                 substitute(paste("(",p == pval,")"), 
                            list(pval=format.pval(w[2], digits=2))),
                 pos=4, cex=par()$cex.axis)
#            legend("topleft", substitute(tau == val, list(val=round(w[1],2) )))
    }
    else if (method=="pearson") {
            text(xshift(xpos), yshift(ypos),
                 substitute(paste(r == val, "\n (p== ", pval, ")"), list(val=round(w[1],2),pval=format.pval(w[2]))), 
            pos=4, cex=par()$cex.lab)
    }
    else if (method=="spearman") {
            text(xshift(xpos), yshift(ypos),
                 substitute(paste(rho == val, " (p== ", pval, ")"), list(val=round(w[1],2),pval=format.pval(w[2]))), 
            pos=4, cex=par()$cex.lab)
    }

}
	
compute_tau <- function(X, indicator, windowsize=length(X)/2, method=c("kendall", "pearson", "spearman"))
## unlike warning_stats, takes indicator as character instead of a function
## assumes X is ts object -- should add to a check(?)
{
	Y <- compute_indicator(X, indicator, windowsize)
	out <- cor.test(time(X)[windowsize:length(X)], Y, method=method)
	c(out$estimate, out$p.value)
}




all_indicators <- function(X, indicators = c("Variance", "Autocorrelation", "Skew", "Kurtosis", "CV"), method=c("kendall", "pearson", "spearman"), ...)
## Calc and plot all the leading indicators in a single frame plot
##		using a simple loop over the plot_indicator fn
## Args 
##		X -- can be a list of n ts objects or a single ts object
##		indicators -- a list of indicators m to calculate
## Returns:
##      Plot m+1 x n matrix of plots, showing data across the first row
##		 and various indicators below.  Columns for different datasets
{
	if(is(X, "list")){
		n <- length(X) # number of datasets
	} else if (is(X, "ts")){
		n <- 1 # number of datasets
		X <- list(X) 
	} else { 
		warning("class of X not recognized")
	}
    data_names <- names(X)
	m <- length(indicators) # number of indicators

## mar is inner margins, in order bottom, left, top, right. 
## oma is outer margins, default to 0 
	par(mfrow=c(m+1,n), oma=c(3,3,2,.2), mar=c(0,1,0,1), ...)
	for(i in 1:n){
		plot(X[[i]], type="l", ylab="Data", xaxt="n", ...)
		mtext(data_names[i],  NORTH<-3, cex=par()$cex.lab, line=1) ## data names on each col
		if(i==1) mtext("Data", WEST<-2, line=3, cex=par()$cex.lab, las=0)  ## "data" y-axis label
	}
## Starts on next row, and goes across datasets
	for(j in 1:m){
		if(j == m){ xaxt <- "s"
		} else {	xaxt <- "n"
		}
		for(i in 1:n){
			plot_indicator(X[[i]], indicators[j], xaxt=xaxt, method=method, xpos=-15, ...) 
			if(i==1) mtext(indicators[j], WEST<-2, line=3, cex=par()$cex.lab, las=0) ## stat name on each row
## i == 2 breaks generality of this plot, making the x-axis appear only on the second column always
			if(j==m & i==2) mtext("Time", SOUTH<-1, line=2, cex=par()$cex.lab) ## x-axis label
		}
	}
}



## a quick labling functions, find percent distance along x and y axis, starting at origin
xshift <- function(xsteps){
	deltax <- (par()$xaxp[2]-par()$xaxp[1])/100
	par()$xaxp[1]+xsteps*deltax
}
yshift <- function(ysteps){
	deltay <- (par()$yaxp[2]-par()$yaxp[1])/100
	par()$yaxp[1]+ysteps*deltay
}





################ DEPRICATED?? ##############

warning_stats <- function(X, indicator, method=c("kendall", "pearson", "spearman")){
	method <- match.arg(method)
	if(is(X,"ts")){
		w <- length(X)/2
		end <- length(X)
		out <- cor.test(time(X)[w:end], indicator(X), method=method)
	} else {
		w <- length(X[,1])/2
		end <- length(X[,1])
		out <- cor.test(X[w:end,1], indicator(X[,2]), method=method)
	}
	c(out$estimate, out$p.value)
}

show_stats <- function(X, indicator, xpos=20, ypos=0, method=c("kendall", "pearson", "spearman")){
		w <- warning_stats(X, indicator)
	text(xshift(xpos), yshift(ypos), paste(method, "coef=", round(w[1],2), "pval=", format.pval(w[2]) ) 
#		 substitute(paste(method, "coef=", val, " (p ", pval, ")"), list(val=round(w[1],2),pval=format.pval(w[2])))
	)
}


