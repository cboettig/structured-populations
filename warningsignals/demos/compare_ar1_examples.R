require(warningsignals)
require(socialR)

# CaCO3 record
caco3 <- read.table("../data/caco3.txt")
# labels, time is in millions of years Before Present, so make negative
caco3 <- data.frame("MYYrs"=-caco3$V1, "CaCO3"=caco3$V2)

## window the data
g_ca <- caco3$MYYrs >= -39 & caco3$MYYrs <= -32  # Data with collapse (for plot)
p_ca <- caco3$MYYrs >= -39 & caco3$MYYrs < -34  # Data used in warning signal
X <- data.frame("time"=caco3$MYYrs[p_ca], "data"=caco3$CaCO3[p_ca])
# Rather annoying to have time backwards and negative, lets reverse this.
X <- data.frame("time"=rev(X[,1] - min(X[,1])), "data"=rev(X[,2]))
## Also annoying that there are replicate values, luckily a quick averaging call will remove them. 
require(limma)
X <-avereps(X, ID=X[,1])
## Time is linearly interpolated for even spacing in the Dakos approach
Y<- approx(X[,1], X[,2], n=482)

## smooth the interopated data.  windowsize default?
smooth <- ksmooth(Y$x, Y$y, bandwidth=.2)
Z <- Y
Z$y <- Y$y-smooth$y

## Transform into a timeseries
start<-Z$x[1]; end<-Z$x[length(Z$x)]
X_ts <- ts(Z$y, start=start, end=end, frequency=length(Z$y)/(end-start))

plt <- function(){
par(mfrow=c(2,1))
w<-round(length(X_ts)/2)
time_window <- time(X_ts)[w:length(X_ts)]
plot(time_window, window_autocorr(X_ts, w), xlim=c(start(X_ts), end(X_ts)), type="l", main="Autocorrelation", xlab="Time", ylab="autocorrelation")
show_stats(X_ts, window_autocorr, ypos=0)
abline(v=time_window[1], lty="dashed")
plot(time_window, window_ar.ols(X_ts, w), xlim=c(start(X_ts), end(X_ts)), type="l", main="Autocorrelation", xlab="Time", ylab="autocorrelation")
show_stats(X_ts, window_ar.ols)
abline(v=time_window[1], lty="dashed")
}

social_plot(plt, tag="stochpop warningsignals CaCO3", file="compare_ar1_examples.png")
