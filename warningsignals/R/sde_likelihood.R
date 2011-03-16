

## Dependencies: NORMT3 in analytic_V, though now depricated.  otherwise all base 

################ Model definitions and functions ###################################
#### Consider running a profiler and moving these into C code for speed ############

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
	if(length(X) == 0){ return(0)}

	# Handle a bunch of cases for possible parameter specification...
	if(is(pars, "numeric")){
	  if(is.null(names(pars))){ 
		  theta1 <- pars[1]; theta2 <- pars[2]; theta3 <- pars[3] 
	  } else if(names(pars)[1] == "theta1"){ 
		  theta1 <- pars["theta1"]; theta2 <- pars["theta2"]; theta3 <- pars["theta3"]
	  } else if(names(pars)[1] == "alpha"){ 
		  theta1 <- pars["alpha"]*pars["theta"]; theta2 <- pars["alpha"]; theta3 <- pars["sigma"]
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
	out <- Inf
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


simulate.OU <- function(pars, t0 = pars$t0, T = pars$T, X0 = pars$X0, N = pars$N){
	if(is.null(t0)) t0 <- 0
	if(is.null(T)) T <- 1
	if(is.null(X0)) X0 <- 1
	if(is.null(N)) N <- 100


	theta = c(theta1=pars$alpha*pars$theta, theta2=pars$alpha, theta3=pars$sigma)
	sde.sim(model="OU", theta= theta, X0=X0, N=N, t0=t0, T=T) 
}

update.OU <- function(pars, X, use_mle=FALSE,method = c("Nelder-Mead",
"BFGS", "CG", "L-BFGS-B", "SANN")){

	method <- match.arg(method)
	if(method == "L-BFGS-B"){ 
		lower <- c(0, -Inf, 0, -Inf) 
	} else { 
		lower = -Inf
	}

	if(!is.null(pars$data)) pars$data <- NULL #needs only parameters 


	if(use_mle){
		OU.mle <- function(alpha, theta, sigma){
			OU.lik(X, list(	alpha=alpha, 
							theta=theta, 
							sigma=sigma))
		}
		fit_results <- mle(OU.mle, 
			start=list(	alpha=pars$alpha,
						theta=pars$theta,
						sigma=pars$sigma), 
			method=method, 
			lower=lower)
	} else {
		X <<- X
		fit_results <- suppressMessages(
				optim(	pars, 
						OU.likfn, 
						method=method, 
						lower=lower, 
						control=list(maxit=2000)))

	}

	if(is(fit_results, "mle")){
		out <- as.list(fit_results@coef)
		out$loglik <- -fit_results@min
		out$T <- time(X)[length(X)]
		out$t0 <- time(X)[1]
		out$X0 <- X[1]
		out$N <- length(X)
		out$data <- X
		class(out) <- class(pars)
	} else {
		out <- as.list(fit_results$par)
		out$loglik <- -fit_results$value
		out$T <- time(X)[length(X)]
		out$t0 <- time(X)[1]
		out$X0 <- X[1]
		out$N <- length(X)
		out$data <- X
		class(out) <- class(pars)
	}
	out
}


########### Early functions, effectively depricated #################


### Early Warning model -- linear change


numeric_V <- function(Dt, pars){
	vint <- function(x){
		exp(-pars$beta*x^2 -2* pars$alpha_0*x)
	}
	int <- integrate(vint, 0, Dt)
	int$value*pars$sigma^2
}

analytic_V <- function(Dt, pars){
	near_zero <- FALSE
	upper_erf <- erf((pars$alpha +Dt*pars$beta)/sqrt(pars$beta) )
	if( upper_erf != 1){ # check numerical stability of erf
		Log_Vx <- log(pars$sigma^2 *  sqrt(pi) ) + pars$alpha_0^2/(pars$beta)  + 
			log( ( upper_erf - erf(pars$alpha_0/(sqrt(pars$beta) ) ))/(2*sqrt(pars$beta) ))
		class(Log_Vx) = "numeric"
		Vx <- exp(Log_Vx)
	} else {
		warning("beta too near zero, using approximation")
		Vx <- pars$sigma^2*(1-exp(-2*pars$alpha*Dt))/(2*pars$alpha)

	}
		Vx
}

# Parametrization dXt = alpha(theta - Xt)dt + sigma*dWt, alpha(t) = beta*t+alpha_0
warning_model <- function(Dt, Xo, t, pars, analytic=FALSE){
	# Assumes pars is a list, as this is the format wanted by mle, optim
    pars$alpha_0 <- pars$alpha_0 + pars$beta*t


	# negative alpha set to near zero
	pars$alpha_0[pars$alpha_0 < 0] <- 1e-10


	int <- pars$beta*Dt^2/2 + pars$alpha_0*Dt
	Ex <- Xo * exp(-int) + pars$theta * (1 - exp(-int) )
	if(analytic)
		Vx <- sapply(1:length(t), function(i){
						params <- pars
						params$alpha_0 <- pars$alpha_0[i]
						analytic_V(Dt, params)	
					})
	else 
		Vx <- sapply(1:length(t), function(i){
						params <- pars
						params$alpha_0 <- pars$alpha_0[i]
						numeric_V(Dt, params)	
					})
    return(list(Ex=Ex,Vx=Vx))
}


dcWarning <- function(x, Dt, x0, pars, t, log = FALSE){
  P <- warning_model(Dt, x0, t, pars)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

rcWarning <- function(n=1, Dt, x0, t, pars){
  P <- warning_model(Dt, x0, t, pars)
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

  if(alpha_0 < 0) return(Inf)
  if(sigma < 0 ) return(Inf)
  pars = list(alpha_0=alpha_0, theta=theta, sigma=sigma, beta=beta)
  -sum( dcWarning(X[2:n], dt, X[1:(n-1)], pars, t=time(X)[1:(n-1)], log=TRUE) )
}

warning.likfn <- function(pars){
			warning.lik(X, pars)
}


simulate.warning <- function(pars, t0 = pars$t0, T = pars$T, X0 = pars$X0, N = pars$N){
	if(is.null(t0)) t0 <- 0
	if(is.null(T)) T <- 1
	if(is.null(X0)) X0 <- 1
	if(is.null(N)) N <- 100

	delta_t <- (T-t0)/N
	Y <- numeric(N)
	Y[1] <- X0
	for(i in 1:(N-1)){
		t <- t0 + (i-1)*delta_t
		Y[i+1] <- rcWarning(1, Dt=delta_t, Y[i], t, pars)
	}
	Y
	ts(Y, start=t0, deltat=delta_t)
}

update.warning <- function(pars, X, use_mle=FALSE,method = c("Nelder-Mead", 
    "BFGS", "CG", "L-BFGS-B", "SANN")){

	method <- match.arg(method)
	if(method == "L-BFGS-B"){ 
		lower <- c(0, -Inf, 0, -Inf) 
	} else { 
		lower = -Inf
	}

	if(!is.null(pars$data)) pars$data <- NULL #needs only parameters 

	if(use_mle){
		warning.mle <- function(alpha_0, theta, sigma, beta){
			warning.lik(X, 
						list(alpha_0=alpha_0,
							theta=theta, 
							sigma=sigma, 
							beta=beta) )
		}
		fit_results <- mle( warning.mle, 
					start=list(	alpha_0=pars$alpha_0, 
								theta=pars$theta, 
								sigma=pars$sigma, 
								beta=pars$beta), 
					method=method, 
					lower=lower)
	} else {
		X <<- X
		fit_results <- optim( pars, 
						warning.likfn, 
						method=method,
						lower=lower,
						control=list(maxit=2000)
						)
	}
	if(is(fit_results, "mle")){
		out <- as.list(fit_results@coef)
		out$loglik <- -fit_results@min
		out$T <- time(X)[length(X)]
		out$t0 <- time(X)[1]
		out$X0 <- X[1]
		out$N <- length(X)
		out$data <- X
		class(out) <- class(pars)
	} else {
		out <- as.list(fit_results$par)
		out$loglik <- -fit_results$value
		out$T <- time(X)[length(X)]
		out$t0 <- time(X)[1]
		out$X0 <- X[1]
		out$N <- length(X)
		out$data <- X
		class(out) <- class(pars)
	}
	out
}


#### ChangePt function library

# pars = c(alpha_1, alpha_2, theta, sigma, t_shift)
simulate.changePt <- function(pars,  t0 = pars$t0, T = pars$T, X0 = pars$X0, N = pars$N){
	if(!is(pars, "list")) message("pars should be list format") 

	if(is.null(t0)) t0 <- 0
	if(is.null(T)) T <- 1
	if(is.null(X0)) X0 <- 1
	if(is.null(N)) N <- 100


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
					deltat=delta_t   ))

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
#	t_shift <- pars$t_shift
#	pars_ou1 <- list(alpha=pars$alpha_1, theta=pars$theta, sigma=pars$sigma)
#	pars_ou2 <- list(alpha=pars$alpha_2, theta=pars$theta, sigma=pars$sigma)

	if(sum(time(X) <= t_shift) == 0){
		X1 <- numeric(0)
	} else {
		X1 <- ts(X[time(X) <= t_shift], start=start(X), deltat=deltat(X))
	}
	if(sum(time(X) > t_shift) == 0){
		 X2 <- numeric(0)
	} else {
		X2 <- ts(X[time(X) > t_shift], start=t_shift, deltat=deltat(X))
	}
	OU.lik(X1, pars_ou1 ) + OU.lik(X2, pars_ou2)
}

changePt.likfn <- function(pars){
	changePt.lik(X, pars) 
}

update.changePt <- function(pars,X,method = c("Nelder-Mead", 
    "BFGS", "CG", "L-BFGS-B", "SANN")){

	method <- match.arg(method)
	if(method == "L-BFGS-B"){ 
		lower <- c(0, 0, -Inf, 0, -Inf) 
	} else { 
		lower = -Inf
	}


	X <<- X
	if(!is.null(pars$data)) pars$data <- NULL #needs only parameters 
	fit_results <- suppressMessages( optim(pars, changePt.likfn, method=method, lower=lower) )

	out <- as.list(fit_results$par)
	out$loglik <- -fit_results$value
	out$T <- time(X)[length(X)]
	out$t0 <- time(X)[1]
	out$X0 <- X[1]
	out$N <- length(X)

	out$data <- X
	class(out) <- class(pars)
	out
}



####################  bootstrapping ##############



bootstrap <- function(model, observed = NULL, reps=4, cpu=2){
	require(snowfall)
	if(! sfIsRunning() ){
		if(cpu<2){
			sfInit()
		} else {
			sfInit(parallel=TRUE, cpu=cpu)
			sfLibrary(stochPop)
			sfExportAll()
		}
	}
	# if data not provided seperately, try and get it from model
	if(is.null(observed)){
		if(!is.null(model$data) ){
			observed <- model$data
		} else { warning("data not found") }
	}


	data <- sfSapply(1:reps, 
		function(i) simulate(model, T=model$T, N=model$N, X0=model$X0) )
	fits <- sfSapply(1:reps, 
		function(i) update(model, data[,i]) )
	
}


plot_bootstrap <- function(object, parameter="all", model=1, ...){
## plot the model specified if given the full likelihood ratio bootstrap, rather than reboostrapping the model
	if(!is(object, "matrix")){
		n_models <- sqrt(length(object$bootstraps))
		j <- model
		object <- object$bootstraps[[(j-1)*n_models+j]]
	}

	if(parameter == "all"){
		n_pars <- dim(object)[1]-6
		par(mfrow=c(1,n_pars) )
		for(i in 1:n_pars){
			plot(density(unlist(object[i,])), xlab=rownames(object)[i], main="", ... )
		}
	} else { 
		i <- pmatch(parameter, rownames(object))
		plot(density(unlist(object[i,])), xlab=rownames(object)[i], main="", ...  )
	}
}


bootstrapLR <- function(model_list, reps=4, cpu=2){
	
	n_models <- length(model_list)
	data <- vector(mode="list", length=n_models)
	fits <- vector(mode="list", length=n_models^2)

	# Calculate the likelihood ratio for the model set
	observed_LR_ratios <- matrix(NA, n_models, n_models)
	for(j in 1:n_models){
		for(k in 1:n_models){
			# 2*(test-null)
			observed_LR_ratios[j,k] <- 2*(model_list[[j]]$loglik - model_list[[k]]$loglik)
		}
	}
	# Initialize parallel processing
	require(snowfall)
	if(cpu<2){
		sfInit()
	} else {
		sfInit(parallel=TRUE, cpu=cpu)
		sfLibrary(stochPop)
		sfExportAll()
	}

	# simulate data from each model
	for(j in 1:n_models){
		data[[j]] <- sfSapply(1:reps, 
					function(i) simulate(model_list[[j]]) 	)
	}
	# fit each model to eack dataset
	for(j in 1:n_models){
		for(k in 1:n_models){
			fits[[(j-1)*n_models+k]] <- sfSapply(1:reps, 
					function(i) update(model_list[[j]], data[[k]][,i]) )
		}
	} 
# consider output of original model values, etc
	out <- list(bootstraps = fits, observed_LR_ratio = observed_LR_ratios)
	class(out) = "LR_bootstraps"
	out
}


# add function to reconstruct the plots 


# plot likelihood ratio of model i vs model j
LRplot <-  function(input, test_index, null_index, ...){
	i<- test_index
	j<- null_index
	object <- input$bootstraps
	n_models <- sqrt(length(object))
	# model i on dataset j
	test <- object[[(i-1)*n_models+j]]
	# model j on dataset j
	null <- object[[(j-1)*n_models+j]]

	test_l <- pmatch("loglik", rownames(test) )
	null_l <- pmatch("loglik", rownames(null) )

	lr <- 2*( unlist( test[test_l,] ) - unlist( null[null_l, ] ) )
	
	obs_lr <- input$observed_LR_ratio[i,j]
	xlim = range(c(range(lr), obs_lr) )
	hist(lr, border='white', col='lightblue', xlim=xlim, ...)
	abline(v=obs_lr , lwd=4, lty=3, col="darkblue")
	text(.98*obs_lr, .5*par()$yaxp[2], paste("p = ", round(sum(lr > obs_lr)/length(lr), digits=4)), cex=1.5)

}


########### Power Test ##############
power_pair <- function(null, test, nboot = 100, cpu = 2, threshold = .95, name_parameter = "beta"){

	## Gotta get templates for the models, do so by fitting some dummy data
	
	## are we in parallel?
	if(cpu>1){ 	
		sfInit(parallel=TRUE, cpu=cpu) 
		sfLibrary(stochPop)
		sfExportAll()
	} else sfInit()

	null_dist <- sfSapply(1:nboot, function(i){
		data <- simulate(null)
		null <- update(null, X=data)
		test <- update(test, X=data)
		-2*(null$loglik - test$loglik) 
	})

	## Actually do the bootstraps of the test model for each alpha


	test_dist <- sfSapply(1:nboot, function(i){
			data <- simulate(test)
			null <- update(null, data)
			test <- update(test, data)
			-2*(null$loglik - test$loglik) 
	})

	## Power calculation
		threshold_tail <- sort(null_dist)[ round(threshold*nboot) ]
		power <- sum(test_dist > threshold_tail)/nboot

	## format the output
	output <- list(	null_dist=null_dist, test_dist=test_dist, power=power, 
					name_parameter=name_parameter, 
					nboot=nboot, threshold=threshold, null=null,test=test)
	class(output) <- "power_pair"
	output
}

plot.power_test <- function(object){
}





power_test <- function(null, test, nboot = 100, cpu = 2, threshold = .95, name_parameter = "beta", values_parameter = seq(0.00001, 100, length=10)){

	## Gotta get templates for the models, do so by fitting some dummy data
	
	## are we in parallel?
	if(cpu>1){ 	
		sfInit(parallel=TRUE, cpu=cpu) 
		sfLibrary(stochPop)
		sfExportAll()
	} else sfInit()

	null_dist <- sfSapply(1:nboot, function(i){
		data <- simulate(null)
		null <- update(null, X=data)
		test <- update(test, X=data)
		-2*(null$loglik - test$loglik) 
	})

	## Actually do the bootstraps of the test model for each alpha
	test_dist <- sfLapply(1:length(values_parameter), function(i){

		test[paste(name_parameter)] <- values_parameter[i]

		sfSapply(1:nboot, function(i){
				data <- simulate(test)
				null <- update(null, data)
				test <- update(test, data)
				-2*(null$loglik - test$loglik) 
		})
	})

	## Power calculation
	power <- sapply(1:length(values_parameter), function(i){
		threshold_tail <- sort(null_dist)[ round(threshold*nboot) ]
		sum(test_dist[[i]] > threshold_tail)/nboot
	})

	## format the output
	output <- list(	null_dist=null_dist, test_dist=test_dist, power=power, 
					name_parameter=name_parameter, values=values_parameter, 
					nboot=nboot, threshold=threshold, null=null,test=test)
	class(output) <- "power_test"
	output
}

plot.power_test <- function(object){
}

