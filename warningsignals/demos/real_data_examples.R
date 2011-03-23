
all_indicators <- function(X, indicators = c("Variance", "Autocorrelation", "Skew", "Kurtosis"))
## Plot all the leading indicators in a single frame plot
## Args 
##		X -- can be a data.frame of n ts objects or a single ts object
##		indicators -- a list of indicators m to calculate
## Returns:
##      Plot m+1 x n matrix of plots, showing data across the first row
##		 and various indicators below.  Columns for different datasets
{
	if(is(X, "data.frame")){
		n <- length(X)
	} else if (is(X, "ts")){
		n <- 1
		X <- data.frame(X)
	} else { 
		warning("class of X not recognized")
	}
	m <- length(indicators)

## mgp is margin of title, axis label, axis line.  3,1,0 is default
#	par(cex.lab=1.7, lwd=2, mgp=c(2,.4,0) )

	par(mfrow=c(m+1,n))

## postitions of the plots 1, 2, 3, 4, 5 in a matrix layout
#	mat <-	rbind(c(1),c(2), c(3), c(4), c(5) )
#	layout(mat, height = c(1.4,1,1,1,1.45))

## mar is margins, in order bottom, left, top, right.  default is 5,4,4,2
	par( mar=c(0,6,4,2), xaxt="n") ## top margin
	for(i in 1:n){
		plot(X[[i]], type="o", ylab="data")
	}

	for(j in 1:m){
		plot_indicator(X[[i]], indicator[j])
	}

	if(m > 1){
		par( mar=c(0,6,0,2), xaxt="n") ## no top or bottom margin
	}


	plot_indicator(X, "Variance", xaxt="n")
	plot_indicator(X, "Autocor", xaxt="n")
	plot_indicator(X, "Skew", xaxt="n")
	# plot_indicator(X, "CV")
	par( mar=c(5,6,0,2) ) ## restore bottom margin
	plot_indicator(X, "Kurtosis")
}













source("load_CaCO3.R")
social_plot(all_indicators(X), tags=tags, file="CaCO3.png")

source("load_deut.R")
## Glaciation I, processed
social_plot(all_indicators(data[[1]]$X_ts), tags=tags, file="DeutI.png")
## Glaciation II, processed
social_plot(all_indicators(data[[2]]$X_ts), tags=tags, file="DeutII.png")
## Glaciation III, processed
social_plot(all_indicators(data[[3]]$X_ts), tags=tags, file="DeutIII.png")
## Glaciation IV, processed
social_plot(all_indicators(data[[4]]$X_ts), tags=tags, file="DeutIV.png")


## Glaciation I, unprocessed
social_plot(all_indicators(data[[1]]$raw_ts), tags=paste(tags, "deut1", "raw"), file="DeutI_raw.png")
## Glaciation II, unprocessed
social_plot(all_indicators(data[[2]]$raw_ts), tags=paste(tags, "deut2", "raw"), file="DeutII_raw.png")
## Glaciation III, unprocessed
social_plot(all_indicators(data[[3]]$raw_ts), tags=paste(tags, "deut3", "raw"), file="DeutIII_raw.png")
## Glaciation IV, unprocessed
social_plot(all_indicators(data[[4]]$raw_ts), tags=paste(tags, "deut4", "raw"), file="DeutIV_raw.png")



##### should be able to pass raw time series in matrix forms, but no luck...
## Glaciation I, unprocessed
social_plot(all_indicators(data[[1]]$X), tags=tags)
## Glaciation II, unprocessed
social_plot(all_indicators(data[[2]]$X), tags=tags)
## Glaciation III, unprocessed
social_plot(all_indicators(data[[3]]$X), tags=tags)
## Glaciation IV, unprocessed
social_plot(all_indicators(data[[4]]$X), tags=tags)


