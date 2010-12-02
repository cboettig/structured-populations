#gaussian_process.R


## Note that setmodel must have form:
## set_a_gauss_model <- function(Dt, Xo, t, pars){
## see the linearized models in likelihood_bifur_models


init_gauss <- function(pars, setmodel, N=100, Xo=1, T=1, t0=0){
# pars are model parameters
# setmodel is a function returning Ex and Vx of the gaussian process
	out <- list(Xo=Xo, t0=t0, T=T, pars=pars, N=N, setmodel=setmodel)
	class(out) <- "gauss"
	out
}
rc.gauss <- function(setmodel, n=1, Dt, x0, t, pars){
	P <- setmodel(Dt, x0, t, pars)
	rnorm(n, mean=P$Ex, sd=sqrt(P$Vx))
}
dc.gauss  <- function(setmodel, x, Dt, x0, t, pars, log = FALSE){
  P <- setmodel(Dt, x0, t, pars)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}
pc.gauss  <- function(setmodel, x, Dt, x0, t, pars, lower.tail = TRUE, log.p = FALSE){ 
  P <- setmodel(Dt, x0, t, pars)
  pnorm(x, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p)
}

lik.gauss <- function(X, pars, setmodel){
	if(length(X) == 0){ return(0)}
    n <- length(X)
    dt <- deltat(X)
	# returns the minus loglik
	out <- -sum( dc.gauss(setmodel, X[2:n], dt, X[1:(n-1)], t=time(X)[1:(n-1)], pars, log=TRUE) )
	out
}
simulateGauss <- function(setmodel, pars, N=100, Xo = 1, T = 1, t0 = 0){
	X <- numeric(N)
	X[1] <- Xo
	delta_t <- (T-t0)/N
	times <- seq(t0, T, by=delta_t)
	for(i in 1:(N-1) ){
		X[i+1] <- rc.gauss(setmodel, 1, Dt=delta_t, x0=X[i], t=times[i+1], pars)
	}
	ts(X, start=t0, deltat=delta_t)
}
updateGauss <- function(setmodel, pars, X, method = c("Nelder-Mead", 
					"BFGS", "CG", "L-BFGS-B", "SANN"), ...){
	method <- match.arg(method)
	likfn <- function(pars) lik.gauss(X, pars, setmodel)
	fit <- optim(pars, likfn, method=method, ...)
	out <- list(pars=fit$par, loglik=-fit$value, T=time(X)[length(X)],
		t0=time(X)[1], Xo <- X[1], X=X, N=length(X), optim_output = fit, setmodel=setmodel)
	class(out) <- "gauss"
	out
}

### For montecarlotest method, need generics: 
simulate.gauss <- function(m) simulateGauss(m$setmodel, pars=m$pars, N=m$N, Xo=m$X[1], T = m$T, t0=m$t0)
update.gauss <- function(m, X, ...) updateGauss(setmodel=m$setmodel, pars=m$pars, X=X, ...)
loglik.gauss <- function(m) m$loglik
getParameters.gauss <- function(m) m$pars


