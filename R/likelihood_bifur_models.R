# set the value of the saddle node bifurcation. 
require(odesolve)
setSN <- function(Dt, Xo, pars){
## Names make it easy to read but slow execution.  For speed, can rewrite in C, see ?lsoda

	Xo <- 1
	pars <- c("r"=1, "theta"=0, "beta"=1)
	Dt = 10
	moments <- function(t,y,p){ 
		yd1 <- pars["r"] - (y[1] - p["theta"])^2 
		yd2 <- -2*(y[1] - p["theta"])*y[2] + p["beta"] 
		list(c(yd1=yd1, yd2=yd2))
	}
	times <- c(0, Dt)
	out <- lsoda( c(Xo, 0), times, moments, pars)
	Ex <- out[2,2] # times are in rows, cols are time, par1, par2
	Vx <- out[2,3]
	return(list(Ex=Ex, Vx=Vx))
}

rcSN <- function(n=1, Dt, x0, pars){
  P <- setOU(Dt, x0, pars)
  rnorm(n, mean=P$Ex, sd = sqrt(P$Vx)) 
}


dcSN <- function(x, Dt, x0, pars, log = FALSE){
  P <- setSN(Dt, x0, pars)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

pcOU <- function(x, Dt, x0, pars, lower.tail = TRUE, log.p = FALSE){ 
  P <- setSN(Dt, x0, pars)
  pnorm(x, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p)
}


# sde likelihood, should take X as an argument and pars as an argument
OU.lik <- function(X, pars){
	if(length(X) == 0){ return(0)}
    n <- length(X)
    dt <- deltat(X)
	out <- -sum( dcSN(X[2:n], dt, X[1:(n-1)], c(theta1, theta2, theta3), log=TRUE) )
	out
}

SN.likfn <- function(pars){
			SN.lik(X, pars)
}



