# likelihood_bifur_models.R

## Dependencies: odesolve
require(odesolve)


## The linearized bifurcation models below use the new gaussian_proccess library for rc, dc, pc, lik, sim, update, etc



###############  Linearized Saddle-Node (LSN) Model ##################
# Will depend explicitly on t
setLSN <- function(Xo, to, t1, pars, R){
	if(any(R(t1,pars) <= 0)){
		Ex <- Xo
		Vx <- rep(Inf,length(Xo))
	} else {
		moments <- function(t,y,p){ 
			sqrtR <- sqrt(R(t,pars)) 
			yd1 <- 2*sqrtR*(sqrtR+pars['theta'] - y[1]) 
			yd2 <- -2*sqrtR*y[2] + p["sigma"]^2*(sqrtR+pars['theta'])
			list(c(yd1=yd1, yd2=yd2))
		}
		jacfn <- function(t,y,p){
			sqrtR <- sqrt(R(t,pars)) 
			c(
			-2*sqrtR, 0,
			0, -2*sqrtR
		)}
## The apply calls needed to work with vector inputs as Xo (whole timeseries)
		times <- matrix(c(to, t1), nrow=length(to))
		out <- lapply(1:length(Xo), function(i){
			lsoda(y=c(xhat=Xo[i], sigma2=0), times=times[i,], func=moments, 
			      parms=pars, jac=jacfn) 
		})
		Ex <- sapply(1:length(Xo), function(i) out[[i]][2,2]) # times are in rows, cols are time, par1, par2
		Vx <- sapply(1:length(Xo), function(i) out[[i]][2,3])
	}
## Handle badly defined parameters by creating very low probability returns, needed particularly
## for the L-BFGS-B bounded method, which can't handle Infs and does rather poorly...
## Worth investigating a better way to handle this.  
## Note errors can apply to some of the timeseries, possibly by estimates of m that
## force system into bifurcation on later parameters (for which terms become undefined)
	if (pars['sigma'] <= 0){
		warning(paste("sigma=",pars["sigma"]))
		Vx <- rep(Inf, length(Xo))
	}
	if (any(Vx < 0, na.rm=TRUE)){
		warning(paste("discard negative Vx,  "))
		Vx[Vx<0] <- Inf
	}
	return(list(Ex=Ex, Vx=Vx))
}



############## Linearized TransCritical (LTC) Model #################
## sample pars for LTC: 
# pars <- c(Ro=1, m=0, theta=1, sigma=1)

## Will depend explicitly on t
setLTC <- function(Xo, to, t1, pars, R){
	moments <- function(t,y,p){
        r <- R(t,pars)
		yd1 <- r*(r -  y[1]/pars['theta']) 
		yd2 <- -2*r*y[2] + (1+r)*p["sigma"]^2 
		list(c(yd1=yd1, yd2=yd2))
	}
	jacfn <- function(t,y,p){
		c(
		-R(t,pars), 0,
		0, -R(t,pars)
	)}
	times <- matrix(c(to, t1), nrow=length(to))
	out <- lapply(1:length(Xo), function(i){
		lsoda(y=c(xhat=Xo[i], sigma2=0), times=times[i,], func=moments, parms=pars, jac=jacfn) 
	})
	Ex <- sapply(1:length(Xo), function(i) out[[i]][2,2]) # times are in rows, cols are time, par1, par2
	Vx <- sapply(1:length(Xo), function(i) out[[i]][2,3])

## Handle badly defined parameters by creating very low probability returns
	if(pars['sigma'] < 0 ) Vx <- rep(Inf, length(Xo)) 
	return(list(Ex=Ex, Vx=Vx))
}

## Just an OU model where alpha depends on time
setOU_timedep <- function(Xo, to, t1, pars, R){
	moments <- function(t,y,p){ 
		yd1 <- R(t,pars)*(pars['theta'] - y[1]) 
		yd2 <- -2*R(t,pars)*y[2] + p["sigma"]^2 ##check?
		list(c(yd1=yd1, yd2=yd2))
	}
	jacfn <- function(t,y,p){
		c(
		-R(t,pars), 0,
		0, -R(t,pars)
	)}
	times <- matrix(c(to, t1), nrow=length(to))
	out <- lapply(1:length(Xo), function(i){
		lsoda(y=c(xhat=Xo[i], sigma2=0), times=times[i,], func=moments, parms=pars, jac=jacfn) 
	})
	Ex <- sapply(1:length(Xo), function(i) out[[i]][2,2]) # times are in rows, cols are time, par1, par2
	Vx <- sapply(1:length(Xo), function(i) out[[i]][2,3])

## Handle badly defined parameters by creating very low probability returns
	if(pars['sigma'] < 0 ) Vx <- rep(Inf, length(Xo)) 
	return(list(Ex=Ex, Vx=Vx))
}



constOU <- function(Xo, to, t1, pars){
	Dt <- t1 - to
	Ex <- pars["theta"]*(1 - exp(-pars["Ro"]*Dt)) + Xo*exp(-pars["Ro"]*Dt) 
	Vx <- 0.5*pars["sigma"]^2 *(1-exp(-2*pars["Ro"]*Dt))/pars["Ro"]
	if(pars['Ro'] < 0 ) Vx <- rep(Inf, length(Xo)) 
	if(pars['sigma'] < 0 ) Vx <- rep(Inf, length(Xo)) 
	return(list(Ex=Ex, Vx=Vx))
}





## the generic bifurcation parameter: here we assume linear model
R <- function(t, pars){pars['Ro'] + pars['m']*t }
timedep_LTC <- function(Xo, to, t1, pars) setLTC(Xo, to, t1, pars, R)
timedep_LSN <- function(Xo, to, t1, pars) setLSN(Xo, to, t1, pars, R)

const_R <- function(t, pars) pars['Ro']
const_LTC <- function(Xo, to, t1, pars) setLTC(Xo, to, t1, pars, const_R)
const_LSN <- function(Xo, to, t1, pars) setLSN(Xo, to, t1, pars, const_R)




#############
#############  Quadratic Saddle Node Bifurcation Model
#############




setSN <- function(Dt, Xo, pars){
	p_tmp <- as.numeric(pars)
	names(p_tmp) <- names(pars)
	pars<-p_tmp
## Names make it easy to read but slow execution.  For speed, can rewrite in C, see ?lsoda
	moments <- function(t,y,p){ 
		yd1 <- p["r"] - (y[1] - p["theta"])^2 
		yd2 <- -2*(y[1] - p["theta"])*y[2] +
			p["beta"]*abs(p["r"] + (y[1] - p["theta"])^2)
		list(c(yd1=yd1, yd2=yd2))
	}
	jacfn <- function(t,y,p){
		c(
		-2*(y[1]-p["theta"]), 0,
		-2*y[2]+p["beta"]*abs(2*(y[1]-p["theta"])), -2*(y[1] - p["theta"])
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
    "BFGS", "CG", "L-BFGS-B", "SANN"), use_mle=FALSE){
	method <- match.arg(method)
	if(method == "L-BFGS-B"){ 
		lower <- c(-Inf, -Inf, 0) 
	} else { 
		lower = -Inf
	}

	if(use_mle){
		SN.mle <- function(r, theta, beta) SN.lik(X, list(r=r, theta=theta, beta=beta))
		start <- list(r=m$pars['r'], theta=m$pars['theta'], beta=m$pars['beta'])
		fit <- mle(SN.mle, start=start, method=method, lower=lower)
	} else {
		X <<- X
		fit <- optim(m$pars, SN.likfn, method=method, lower=lower)
	}
	if(is(fit, "mle")){
		out <- fit
	} else {
		out <- list(pars=fit$par, loglik=-fit$value, T=time(X)[length(X)],
		t0=time(X)[1], Xo <- X[1], data=X, N=length(X), optim_output = fit)
	}
	class(out) <- m$model
	out
}







#############
#############  Misc useful functions (plotting, etc)
#############


## stable root of the SN model
stableroot <- function(m){
	p <- m$p
	xhat <- sqrt(p['r']) + p['theta']
	names(xhat) <- 'xhat'
	xhat
}



## Specify the plotting function for the quadratic curves,
plotcurves <- function(m, out, pars, X, oucurve=NULL){ # m is the initial conditions model (input to update.SN), out the fit model (output of update.SN), pars the parameters used for the simulated data, X
	xmin <- 0; xmax <- 5+stableroot(out)	
curve( -(x-pars['theta'])^2+pars['r'], xmin, xmax, ylim=c(-2, pars['r']+1), lwd=3, main="true vs estimated model", col="darkgray")
	curve( -(x-out$pars["theta"])^2+out$pars["r"], xmin, xmax, add=T, col="red", lty=2, lwd=3)
	if(!is.null(oucurve)){
		curve( -(x-oucurve$pars["theta"])^2+oucurve$pars["r"], xmin, xmax, add=T, col="green", lty=2, lwd=3)
	}
#text(1,4, paste("N = ", m$N, " T = ", m$T),pos=4)
	text(1, 3, paste("est: ", "r = ", as.character(round(out$pars[1],2)), "theta = ", as.character(round(out$pars[2],2)), "beta = ", as.character(round(out$pars[3],2))), pos=4)
	text(1, 2, paste("init: ", "r = ", as.character(round(m$pars[1],2)), "theta = ", as.character(round(m$pars[2],2)), "beta = ", as.character(round(m$pars[3],2))), pos=4)
	text(1, 1, paste("true: ", "r = ", as.character(round(pars[1],2)), "theta = ", as.character(round(pars[2],2)), "beta = ", as.character(round(pars[3],2))), pos=4)
	abline(h=0, lty=2)
	par(new=TRUE)
	plot(X@.Data, time(X),,type="l",col=rgb(0,0,1,.4),xlim=c(xmin,xmax), xaxt="n",yaxt="n",xlab="",ylab="")
	axis(4)
	mtext("data",side=4,line=3)
}






