#indicator_vs_likelihood.R
tags <- "warningsignals stochpop"
nboot <- 16
require(socialR)
require(warningsignals)
sfInit(parallel=TRUE, cpu=16)


sfLibrary(warningsignals)

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




########################### Begin actual analysis ######################## 




pars <- c(Ro=5.0, m= -.04, theta=100, sigma=1)
const_pars <- c(Ro=5.0, theta=100, sigma=1)
sfExportAll()

## Look at the distribution of Taus
test_tau_dist_var <- sfSapply(1:nboot, function(i){
	X <- simulateGauss(timedep_LTC, pars, N=500, T=100, Xo=100)
	warning_stats(X, window_var)
})
null_tau_dist_var <- sfSapply(1:nboot, function(i){
	Y <- simulateGauss(const_LTC, const_pars, N=500, T=100, Xo=100)
	warning_stats(Y, window_var)
})
save(list=ls(), file="indicator_vs_likelihood.Rdat")
social_plot(plt_tau(test_tau_dist_var, null_tau_dist_var, "Variance"), file="taudist_var.png", tags="warningsignals stochpop tau variance")

## Look at the distribution of Taus on autocorrelation
test_tau_dist_acor <- sfSapply(1:nboot, function(i){
	X <- simulateGauss(timedep_LTC, pars, N=500, T=100, Xo=100)
	warning_stats(X, window_autocorr)
})
null_tau_dist_acor <- sfSapply(1:nboot, function(i){
	Y <- simulateGauss(const_LTC, const_pars, N=500, T=100, Xo=100)
	warning_stats(Y, window_autocorr)
})
save(list=ls(), file="indicator_vs_likelihood.Rdat")
social_plot(plt_tau(test_tau_dist_acor, null_tau_dist_acor, "Autocorrelation"), file="taudist_autcorr.png", tags="warningsignals stochpop tau autocorr")


## Simulate some sample data under slow linear change and under no change.  
warning <- simulateGauss(timedep_LTC, pars, N=500, T=100, Xo=100)
no_warning <- simulateGauss(const_LTC, const_pars, N=500, T=100, Xo=100)


save(list=ls(), file="indicator_vs_likelihood.Rdat")
social_plot(plt_data(warning, no_warning), file="indicators.png", tags=tags)


# Likelihood Fits to each data-set and their relative model comparison
timedep <- updateGauss(timedep_LTC, pars, warning, control=list(maxit=1000))
const <- updateGauss(const_LTC, const_pars, warning, control=list(maxit=1000))
llik_warning_fit <- 2*(loglik(timedep)-loglik(const))
## Ideally we want these to be exact, not estimated. since starting with optimal values, they shouldn't have drifted away on expectation, but MLE is biased estimator (overestimates likelihood)
## This is not of course possible with real data, which must just use the MLEs here.  
timedep$par <- pars
const$par <- const_pars
llik_warning <- 2*(loglik(timedep)-loglik(const))
print(paste("with warning signal:", "LR of fit models ", llik_warning_fit, "LR of true models " llik_warning))


## MONTECARLO Non-parametric bootstrap using the exact values.  as noted above, in real data we would use the MLEs, which has the point-estimate problem.  
out <- montecarlotest(const, timedep, cpu=16, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="indicator_vs_likelihood.Rdat")
social_plot(plot(out), file="indicator_vs_likelihood_mc.png", tag="warningsignal stochpop LTC")



### just compare for random example that LR should be quite small when estimated from the null data
timedep_no <- updateGauss(timedep_LTC, pars, no_warning, control=list(maxit=1000))
const_no <- updateGauss(const_LTC, const_pars, no_warning, control=list(maxit=1000))
llik_nowarning_fit <- 2*(loglik(timedep_no)-loglik(const_no))

## using the exact models (m=0), (const R) should result in an log LR of 0
timedep_no$par <- pars
timedep_no$par['m'] = 0
const_no$par <- const_pars
llik_nowarning <- 2*(loglik(timedep)-loglik(const))
print(paste("without warning signal:", "LR of fit models ", llik_nowarning_fit, "LR of true models " llik_nowarning))

save(list=ls(), file="indicator_vs_likelihood.Rdat")


