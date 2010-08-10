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




# sde likelihood fits exploration
OU.lik <- function(theta1, theta2, theta3){
  if( theta2 < 0 | theta3 < 0){
    out <- 1e12 
  } else {
    n <- length(X)
    dt <- deltat(X)
   out <- -sum( dcOU(X[2:n], dt, X[1:(n-1)], c(theta1, theta2, theta3), log=TRUE) )
  }
  out
}



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



warning.lik <- function(alpha_0, theta, sigma, beta){
  n <- length(X)
  dt <- deltat(X)
  pars = list(alpha_0=alpha_0, theta=theta, sigma=sigma, beta=beta)
  -sum( dcWarning(X[2:n], dt, X[1:(n-1)], pars, log=TRUE) )
}



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
warning.sim <- function(t0 = 0, T = 1, X0 = 1, N = 100, pars ){
	delta_t <- (T-t0)/N
	Y <- numeric(N)
	for(i in 1:(N-1)){
		Y[i+1] <- rcWarning(1, Dt=delta_t, Y[i], pars)
	}
	Y
}



