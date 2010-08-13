# Parametrization dXt = (O1 - O2*Xt)dt + O3*dWt
# Some useful functions from the sde library that are not exported:
checkOU <- function(theta){
  if(theta[2]<=0) stop("the process is not stationary")
  if(theta[3]<=0) stop("variance must be positive")
}

setOU <- function(Dt, x0, theta){
  Ex <- theta[1]/theta[2]+(x0-theta[1]/theta[2])*exp(-theta[2]*Dt)
  Vx <- theta[3]^2*(1-exp(-2*theta[2]*Dt))/(2*theta[2])
  return(list(Ex=Ex,Vx=Vx))
}

rcOU <- function(n=1, Dt, x0, theta){
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  rnorm(n, mean=P$Ex, sd = sqrt(P$Vx)) 
}


dcOU <- function(x, Dt, x0, theta, log = FALSE){
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

pcOU <- function(x, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE){ 
  checkOU(theta)
  P <- setOU(Dt, x0, theta)
  pnorm(x, mean=P$Ex, sd=sqrt(P$Vx),
	lower.tail = lower.tail, log.p = log.p)
}




# sde likelihood, should take X as an argument and pars as an argument
OU.lik <- function(X, pars){
	# Handle a bunch of cases for possible parameter specification...
	if(is(pars, "numeric")){
	  if(is.null(names(pars))){ theta1 <- pars[1]; theta2 <- pars[2]; theta3 <- pars[3] 
	  } else if(names(pars)[1] == "theta1"){ theta1 <- pars["theta1"]; theta2 <- pars["theta2"]; theta3 <- pars["theta3"]
	  } else if(names(pars)[1] == "alpha"){ theta1 <- pars["alpha"]*pars["theta"]; theta2 <- pars["alpha"]; theta3 <- pars["sigma"]
	  }
	} else if(is(pars, "list")){
	  if(is.null(names(pars))){ theta1 <- pars[[1]]; theta2 <- pars[[2]]; theta3 <- pars[[3]] 
	  } else if(names(pars)[1] == "theta1"){ theta1 <- pars$theta1; theta2 <- pars$theta2; theta3 <- pars$theta3 
	  } else if(names(pars)[1] == "alpha"){ theta1 <- pars$alpha*pars$theta; theta2 <- pars$alpha; theta3 <- pars$sigma 
	  }
	}
  # throw large -loglik in case parameters are negative
  if( theta2 <= 0 | theta3 <= 0){
	message("certain parameters cannot be negative")
	out <- 1e12
  # Actually compute the minus log likelihood
  } else {
    n <- length(X)
    dt <- deltat(X)
   out <- -sum( dcOU(X[2:n], dt, X[1:(n-1)], c(theta1, theta2, theta3), log=TRUE) )
  }
  out
}

OU.likfn <- function(pars){
			OU.lik(X, pars)
}


OU.sim <- function(t0 = 0, T = 1, X0 = 1, N = 100, pars ){
	theta = c(theta1=pars$alpha*pars$theta, theta2=pars$alpha, theta3=pars$sigma)
	sde.sim(model="OU", theta= theta, X0=X0, N=N, t0=t0, T=T) 
}

OU.fitML <- function(X, pars, use_mle=FALSE){
	if(use_mle){
		OU.mle <- function(alpha, theta, sigma){
			OU.lik(X, list(	alpha=alpha, 
							theta=theta, 
							sigma=sigma))
		}
		out <- mle(OU.mle, 
			start=list(	alpha=pars$alpha,
						theta=pars$theta,
						sigma=pars$sigma), 
			method="L-BFGS-B", 
			lower=c(0,-Inf,0))
	} else {
		X <<- X
		out <- suppressMessages(
				optim(	pars, 
						OU.likfn, 
						method="L-BFGS-B", 
						lower=c(0,-Inf,0)))

	}
	out
}

### Early Warning model -- linear change


numeric_V <- function(Dt, pars){
	vint <- function(x){
		exp(-pars$beta*x^2 -2* pars$alpha_0*x)
	}
	int <- integrate(vint, 0, Dt)
	int$value*pars$sigma^2
}

analytic_V <- function(Dt, pars){
	upper_erf <- erf((pars$alpha +Dt*pars$beta)/sqrt(pars$beta) )
	if( upper_erf != 1){ # check numerical stability of erf
		Log_Vx <- log(pars$sigma^2 *  sqrt(pi) ) + pars$alpha_0^2/(pars$beta)  + log( ( upper_erf - erf(pars$alpha_0/(sqrt(pars$beta) ) ))/(2*sqrt(pars$beta) ))
		class(Log_Vx) = "numeric"
		Vx <- exp(Log_Vx)
	} else {
		warning("beta near zero, using approximation")
		Vx <- pars$sigma^2*(1-exp(-2*pars$alpha*Dt))/(2*pars$alpha)

	}
		Vx
}

# Parametrization dXt = alpha(theta - Xt)dt + sigma*dWt, alpha(t) = beta*t+alpha_0
warning_model <- function(Dt, Xo, pars, analytic=FALSE){
	# Assumes pars is a list, as this is the format wanted by mle, optim
	int <- pars$beta*Dt^2/2 + pars$alpha_0*Dt
	Ex <- Xo * exp(-int) + pars$theta * (1 - exp(-int) )
	if(analytic)
		Vx <- analytic_V(Dt, pars)
	else 
		Vx <- numeric_V(Dt, pars)	
    return(list(Ex=Ex,Vx=Vx))
}


dcWarning <- function(x, Dt, x0, pars, log = FALSE){
  P <- warning_model(Dt, x0, pars)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

rcWarning <- function(n=1, Dt, x0, pars){
  P <- warning_model(Dt, x0, pars)
  rnorm(n, mean=P$Ex, sd=sqrt(P$Vx))
}



warning.lik <- function(X, pars){
  if(is(pars, "numeric")){
	  alpha_0 = pars[1]; theta = pars[2]; sigma=pars[3]; beta = pars[4]
  } else if(is(pars, "list") ){
	  alpha_0 = pars[[1]]; theta = pars[[2]]; sigma=pars[[3]]; beta = pars[[4]]
  }
  n <- length(X)
  dt <- deltat(X)

  if(alpha_0 < 0) return(1e12)
  if(sigma < 0 ) return(1e12)
  pars = list(alpha_0=alpha_0, theta=theta, sigma=sigma, beta=beta)
  -sum( dcWarning(X[2:n], dt, X[1:(n-1)], pars, log=TRUE) )
}

warning.likfn <- function(pars){
			warning.lik(X, pars)
}


warning.sim <- function(t0 = 0, T = 1, X0 = 1, N = 100, pars ){
	delta_t <- (T-t0)/N
	Y <- numeric(N)
	for(i in 1:(N-1)){
		Y[i+1] <- rcWarning(1, Dt=delta_t, Y[i], pars)
	}
	Y
}

warning.fitML <- function(X, pars, use_mle=FALSE){

	if(use_mle){
		warning.mle <- function(alpha_0, theta, sigma, beta){
			warning.lik(X, 
						list(alpha_0=alpha_0,
							theta=theta, 
							sigma=sigma, 
							beta=beta) )
		}
		out <- mle( warning.mle, 
					start=list(	alpha_0=pars$alpha_0, 
								theta=pars$theta, 
								sigma=pars$sigma, 
								beta=pars$beta), 
					method="L-BFGS-B", 
					lower=c(0, -Inf, 0, -Inf))
	} else {
		X <<- X
		out <- optim(	pars, 
						warning.likfn, 
						method="L"
						#lower=c(0,-Inf,0,-Inf)
						)
	}
	out
}


#### ChangePt function library

# pars = c(alpha_1, alpha_2, theta, sigma, t_shift)
changePt.sim <- function(t0= 0, T = 1, X0 = 1, N = 100, pars ){
	if(!is(pars, "list")) error("pars should be list format") 

	delta_t <- (T-t0)/N
	N1 <- floor( (pars$t_shift - t0)/delta_t)
	if(pars$t_shift < T & pars$t_shift > t0){
		part1 <- sde.sim(	model="OU", 
							theta= c(	pars$theta*pars$alpha_1,
										pars$alpha_1,
										pars$sigma), 
							X0=X0, 
							t0=t0,
							N = N1,
							T=pars$t_shift)
		# note that returns start condition, which is the same as the previous interval end, so interval overlaps by this point. 
		part2 <- sde.sim(	model="OU", 
							theta= c(	pars$theta*pars$alpha_2,
										pars$alpha_2,
										pars$sigma), 
							X0=part1[length(part1)], 
							N = N-N1,
							t0=pars$t_shift, 
							T=T)
		return( ts(
					c(part1, part2[2:length(part2)]), 
					start=t0, 
					deltat=deltat(part1)   ))

	## Evaluate and warn if t_shift is outside time interval
	} else if(pars$t_shift > T) {
		message("time of shift occurs after end of time interval requested, returning all regime 1")
		return( sde.sim(	model="OU", 
							theta= c(pars$theta*pars$alpha_1,
									pars$alpha_1,
									pars$sigma), 
							X0=X0, 
							N=N, 
							T=T  ))
	} else {
		message("time of shift occurs before start of time interval requested, returning all regime 2")
		return(  sde.sim(	model="OU", 
							theta= c(pars$theta*pars$alpha_2,
									pars$alpha_2,
									pars$sigma), 
							X0=X0, 
							N=N, 
							T=T	))
	}
}



changePt.lik <- function(X, pars){

	if(is(pars, "list")){ 
		tmp <- unlist(pars)
		names(tmp) <- names(pars)
		pars <- tmp
	}
	t_shift <- pars[5]
	pars_ou1 <- list(alpha=pars[1], theta=pars[3], sigma=pars[4])
	pars_ou2 <- list(alpha=pars[2], theta=pars[3], sigma=pars[4])
	OU.lik( X[time(X) <= t_shift], pars_ou1 ) + OU.lik( X[time(X)>t_shift], pars_ou2)
}

changePt.likfn <- function(pars){
	changePt.lik(X, pars) 
}

changePt.fitML <- function(X,pars){
	X <<- X
	suppressMessages( optim(pars, changePt.likfn, method="L-BFGS-B", lower=c(0,0,-Inf,0,0)) )
}













############## Depricated ######################

# time shift
# pars <- list(alpha_1 = , alpha_2 = , theta = , sigma = , t_shift = )
TwoRates <- function(Dt, Xo, pars){

	gamma <- pars$alpha_1*min(Dt, pars$t_shift) + pars$alpha_2*max(Dt-pars$t_shift,0)

	omega <- pars$theta*( max(0,exp(pars$alpha_2*Dt) - exp(pars$alpha_2*pars$t_shift)))  + 
			 pars$theta*(exp(pars$alpha_1*min(Dt, pars$t_shift) ) - 1)


	if(omega >= Inf){
		message("assuming equilibrium...")
		return(list(Ex=pars$theta, Vx=pars$sigma/pars$alpha_2) )
	}

	Ex <- exp(-gamma)*(Xo + omega) 

	Vx <- exp(-2*gamma) * ( (pars$sigma)^2/(2*pars$alpha_2) * ( max(0, exp( 2*pars$alpha_2*Dt ) - exp(2*pars$alpha_2*pars$t_shift) ) ) +
		  (pars$sigma)^2/(2*pars$alpha_1) * (exp(2*pars$alpha_1*min(Dt, pars$t_shift)) -1)) 
				
    return(list(Ex=Ex,Vx=Vx))
}

dcTwoRates <- function(x, Dt, x0, pars, log = FALSE){
	P <- TwoRates(Dt, x0, pars)
#	print(P$Vx)
#	print(unlist(pars))
	dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

TwoRates.lik <- function(alpha_1, alpha_2, theta, sigma, t_shift){
	n<- length(X)
	dt <- deltat(X)
	pars <- list(alpha_1 = alpha_1, alpha_2 = alpha_2, theta = theta, sigma = sigma, t_shift=t_shift)
	-sum( dcTwoRates(X[2:n], dt, X[1:(n-1)], pars, log=TRUE) )

}

rcTwoRates <- function(n=1, Dt, x0, pars){
	P <- TwoRates(Dt, x0, pars)
	rnorm(n, mean=P$Ex, sd = sqrt(P$Vx)) 
}

TwoRates.sim <- function(t0 = 0, T = 1, X0 = 1, N = 100, pars ){
	delta_t <- (T-t0)/N
	Y <- numeric(N)
	for(i in 1:(N-1)){
		Y[i+1] <- rcTwoRates(1, Dt=delta_t, Y[i], pars)
	}
	Y
}

