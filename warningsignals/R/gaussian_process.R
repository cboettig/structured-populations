##gaussian_process.R

## setmodel is a function handle returning Ex and Vx of the gaussian process
## the possible setmodel functions are defined in likelihood_bifur_models.R
## Note that setmodel must have form:
## set_a_gauss_model <- function(Xo, to, t1, pars){
## see the linearized models in likelihood_bifur_models
##   OLD FORM depended on fixed Dt:  set_a_gauss_model <- function(Dt, Xo, t, pars){


init_gauss <- function(pars, setmodel, N=100, Xo=1, T=1, t0=0){
# pars are model parameters
	out <- list(Xo=Xo, t0=t0, T=T, pars=pars, N=N, setmodel=setmodel, k=length(pars))
	class(out) <- "gauss"
	out
}
rc.gauss <- function(setmodel, n=1, x0, to, t1, pars){
	P <- setmodel(x0, to, t1, pars)
	rnorm(n, mean=P$Ex, sd=sqrt(P$Vx))
}
dc.gauss  <- function(setmodel, x, x0, to, t1, pars, log = FALSE){
	  P <- setmodel(x0, to, t1, pars)
	  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}
pc.gauss  <- function(setmodel, x, x0, to, t1, pars, lower.tail = TRUE, log.p = FALSE){ 
  P <- setmodel(x0, to, t1, pars)
  pnorm(x, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p)
}



### Time can be specified directly in X as ts object or as the first column of X, and data in the second column
lik.gauss <- function(X, pars, setmodel, log=TRUE){
	if(is(setmodel, "character")) return(heteroOU(X,pars))
## Returns minus log likelihood
	if(length(X) == 0){ return(0)}
	if(is(X,"ts")){
	    n <- length(X)
		dt <- deltat(X)
		# returns the minus loglik
		to <- time(X)[1:(n-1)]
		out <- -sum( dc.gauss(setmodel, X[2:n], X[1:(n-1)], to=to, t1=to+dt, pars, log=log) )
	} else {
#		if(length(X) != 2) stop("Must specify X as a ts object or matrix with time in col 1 and data in col 2 ")
	    n <- length(X[,1])
		times <- X[,1]
		Y <- X[,2] ## data values in second column
		
		out <- -sum( dc.gauss(setmodel, Y[2:n], Y[1:(n-1)], to=times[1:(n-1)], t1=times[2:n], pars, log=log) )
	}
	out
}
simulateGauss <- function(setmodel, pars, N=100, Xo = 1, T = 1, t0 = 0, times=NULL){
	## Return a timeseries if times not specified explicitly
	if(is.null(times)){
		X <- numeric(N)
		X[1] <- Xo
	delta_t <- (T-t0)/N
		times <- seq(t0, T, by=delta_t)
		for(i in 1:(N-1) ){
			X[i+1] <- rc.gauss(setmodel, 1, x0=X[i], to=times[i+1],t1=times[i+2], pars)
		}
		out <- ts(X, start=t0, deltat=delta_t)
	} else {
		N <- length(times)
		X <- numeric(N-1)
		X[1] <- Xo
		for(i in 1:(N-2) ){
			X[i+1] <- rc.gauss(setmodel, 1, x0=X[i], to=times[i+1],t1=times[i+2], pars)
		}
		out <- data.frame(times=times[1:(N-1)], X=X)
	}
	out
}

### Time can be specified directly in X as ts object or as the first column of X, and data in the second column.  
updateGauss <- function(setmodel, pars, X, method = c("Nelder-Mead", 
					"BFGS", "CG", "L-BFGS-B", "SANN"), ...){
	method <- match.arg(method)
	likfn <- function(pars) lik.gauss(X, pars, setmodel) #optim routine needs something that is just a function of the parameters to be searched (not a fn of data X)
	fit <- optim(pars, likfn, method=method, ...)

	if(is(X, "ts")){
		out <- list(pars=fit$par, loglik=-fit$value, T=time(X)[length(X)],
			t0=time(X)[1], Xo <- X[1], X=X, N=length(X), optim_output = fit, setmodel=setmodel, k=length(pars), times=NULL)
	} else {
		out <- list(pars=fit$par, loglik=-fit$value, T=X[length(X[,1]), 1],
			t0=X[1,1], Xo <- X[1,2], X=X, N=length(X), optim_output = fit, setmodel=setmodel, k=length(pars), times=X[,1])
	}
	class(out) <- "gauss"
	out
}

### For montecarlotest method, need generics: 
simulate.gauss <- function(m){
	## X is a ts, use method
	if(is.null(m$times)) simulateGauss(m$setmodel, pars=m$pars, N=m$N, Xo=m$X[1], T = m$T, t0=m$t0)
	else simulateGauss(m$setmodel, pars=m$pars, N=m$N, Xo=m$X[1], T = m$T, t0=m$t0, times=m$times)
}
update.gauss <- function(m, X, method = c("Nelder-Mead", 
					"BFGS", "CG", "L-BFGS-B", "SANN"), ...){
	updateGauss(setmodel=m$setmodel, 
		pars=m$pars, X=X, method=method, ...)
}
loglik.gauss <- function(m) m$loglik
getParameters.gauss <- function(m) c(m$pars, convergence=as.numeric(m$optim_output$convergence))
# update the loglik (calculation rather than lookup)
loglik_calc.gauss <- function(m) -lik.gauss(m$X, m$pars, m$setmodel)

