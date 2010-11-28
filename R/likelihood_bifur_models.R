# set the value of the saddle node bifurcation. 
require(odesolve)
setSN <- function(Dt, Xo, pars){
	p_tmp <- as.numeric(pars)
	names(p_tmp) <- names(pars)
	pars<-p_tmp
## Names make it easy to read but slow execution.  For speed, can rewrite in C, see ?lsoda
	moments <- function(t,y,p){ 
		yd1 <- p["r"] - (y[1] - p["theta"])^2 
		yd2 <- -2*(y[1] - p["theta"])*y[2] +
			p["beta"]*abs(p["r"] - (y[1] - p["theta"])^2)
		list(c(yd1=yd1, yd2=yd2))
	}
	jacfn <- function(t,y,p){
		c(
		-2*(y[1]-p["theta"]), 0,
		-2*y[2]+p["beta"]*abs(-2*(y[1]-p["theta"])), -2*(y[1] - p["theta"])
	)
	}
	times <- c(0, Dt)
	out <- lapply(Xo, function(x0){lsoda(y=c(xhat=x0, sigma2=0), times=times, func=moments, parms=pars, jac=jacfn) 
	})
	Ex <- sapply(1:length(Xo), function(i) out[[i]][2,2]) # times are in rows, cols are time, par1, par2
	Vx <- sapply(1:length(Xo), function(i) out[[i]][2,3])
	return(list(Ex=Ex, Vx=Vx))
}

rcSN <- function(n=1, Dt, x0, pars){
  P <- setSN(Dt, x0, pars)
  rnorm(n, mean=P$Ex, sd = sqrt(P$Vx)) 
}


dcSN <- function(x, Dt, x0, pars, log = FALSE){
  P <- setSN(Dt, x0, pars)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

pcSN <- function(x, Dt, x0, pars, lower.tail = TRUE, log.p = FALSE){ 
  P <- setSN(Dt, x0, pars)
  pnorm(x, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p)
}


# sde likelihood, should take X as an argument and pars as an argument
SN.lik <- function(X, pars){
	if(length(X) == 0){ return(0)}
    n <- length(X)
    dt <- deltat(X)
	# returns the minus loglik
	out <- -sum( dcSN(X[2:n], dt, X[1:(n-1)], pars, log=TRUE) )
	out
}

SN.likfn <- function(pars){
			SN.lik(X, pars)
}


simulate.SN <- function(m){
	X <- numeric(m$N)
	X[1] <- m$Xo
	delta_t <- (m$T-m$t0)/m$N
	for(i in 1:(m$N-1) ){
#		t <- m$t0 +(i-1)*delta_t ## absolute time doesn't matter on this fn
		X[i+1] <- rcSN(1, Dt=delta_t, X[i], m$pars)
	}
	ts(X, start=m$t0, deltat=delta_t)
}

update.SN <- function(m, X, method = c("Nelder-Mead", 
    "BFGS", "CG", "L-BFGS-B", "SANN")){
	method <- match.arg(method)
	X <<- X
	fit <- suppressMessages(
		optim(m$pars, SN.likfn, method=method)
	)
	out <- list(pars=fit$par, loglik=-fit$value, T=time(X)[length(X)],
	t0=time(X)[1], Xo <- X[1], data=X, N=length(X))
	class(out) <- c(m$model, "sdemodel")
	out
}




## Generic functions
simulate.sdemodel <- function (m) {
	if(m$model=="SN"){
		sim <- simulate.SN(m)
	}
	sim
}

init_sdemodel <- function(pars, model, N=100, Xo=1, T=1, t0=0){
	out <- list(Xo=Xo, t0=t0, T=T, pars=pars, N=N, model=model)
	class(out) <- c(model, "sdemodel")
	out
}
