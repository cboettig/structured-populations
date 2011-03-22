## Dependencies: limma for avereps function


dakos_data_processing <- function(X, detrend=TRUE, interpolate=TRUE, n=NULL, bw=NULL){
## Description: mimics the data processing of Dakos et al 2008 PNAS.  
##				Linearly interpolaes evenly spaced data & subtracts a Gaussian smoothed 
##              kernel, and turns into a timeseries (ts) object if interpolated to even spacing
## Args:
##		X -- a matrix or data.frame with times in first column and observations in second.  
##           may have replicate observations at a given time, which will be averaged out
##		detrend -- logical.  If gaussian kernel smoothing should be done
##		interpolate -- logical.  If interpolate is false, returns only a matrix (ts requires even spacing)
##		n --The number of points to be used in interpolation.  If not given, will match length of original data
## Returns:
##    X_ts -- a timeseries object that has been detrended
##		Z  -- the detrended data as a matrix 
##		Y -- the interpolated data matrix
##		X -- the original data with replicates averaged out (matrix)
##		smooth -- the smoothed "trend" that has been removed
##		raw_ts -- X data rendered as a timeseries object

## Also annoying that there are replicate values, luckily a quick averaging call will remove them. 
	X <-avereps(X, ID=X[,1])
## Time is linearly interpolated for even spacing in the Dakos approach
	if(interpolate){ 
		if(is.null(n)) n <- length(X[,1]) 
		Y<- approx(X[,1], X[,2], n=n) 
	} else { Y <- X }
## smooth the interopated data.  windowsize default?
	if(detrend){
		if(is.null(bw)) bw <- bw.nrd(Y$x)
		smooth <- ksmooth(Y$x, Y$y, bandwidth=bw, n.points=length(Y$x), x.points=Y$x)
		Z <- Y
		Z$y <- Y$y-smooth$y
	} else { Z <- Y }
## Transform into a timeseries
	start<-Z$x[1]; end<-Z$x[length(Z$x)]
	X_ts <- ts(Z$y, start=start, end=end, frequency=length(Z$y)/(end-start))
## make an unprocessed timeseries
	start<-X[1,1]; end<-X[length(X[,1]),1]
	raw_ts <- ts(X[,2], start=start, end=end, frequency=length(X[,2])/(end-start))

	out <- list(X_ts=X_ts, Z=Z, Y=Y, X=X, smooth=smooth, raw_ts= raw_ts)
	class(out) <- "dakos"
	out
}


## Make a Dakos-style plot
plot.dakos <- function(dat){
	X_ts <- dat$X_ts; X <- dat$X; Y <- dat$Y; Z <- dat$Z; smooth <- dat$smooth
	par(mfrow=c(2,1))
	plot(X, type="o", col="blue", pch='.', cex=3) # Raw data
	points(Y, col="red", pch='.', cex=3) # Interpolated data
	lines(smooth, col="darkgray", lwd=3) # smoothing function
	w<-round(length(X_ts)/2)
	time_window <- time(X_ts)[w:length(X_ts)]
	plot(time_window, window_autocorr(X_ts, w), xlim=c(start(X_ts), end(X_ts)), type="l", main="Autocorrelation", xlab="Time", ylab="autocorrelation")
	abline(v=time_window[1], lty="dashed")
	show_stats(X_ts, window_autocorr)

## Should have ability label axis in original MYrs BP units
}
