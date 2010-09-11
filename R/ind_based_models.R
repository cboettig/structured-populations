# @file: ind_based_models.R
# @author: Carl Boettiger, <cboettig@gmail.com>
# @section DESCRIPTION wrapper to the C code containing the gillespie simulation for individual-based models.  


#	Pars = {n, e, a, K, h, i, Da,Dt} 
# inits[8] = {572, .5, 160, 1000, 200, 0, 1, 100};

saddle_node_ibm <- function(
	# Xo is initial population size
	# e is the natural per-capita death-rate
	# a is the environmental toxin level
	# K scales the birth rate (hence the equilibrium size)
	# h is the half-max growth rate
	# i is a place-holder for an internal counter, not real parameter
	# Da is the rate of environmental degradation
	# Dt is the time at which environmental degradation begins
	pars=c("Xo" = 570, "e" = .5, "a" = 160, "K" = 1000, "h" = 200, "i" = 0, "Da" = 1, "Dt" = 100),
	times = seq(0,150,length=50), 
	reps=1)
{
	samples <- length(times)
	N <- reps*samples
	maxtime <- max(times)

	o <- .C("saddle_node_direct", double(N), as.double(pars), as.integer(samples), as.integer(reps), as.double(maxtime) )

	x1 = matrix(o[[1]], samples, reps)
	m1 <- sapply(1:samples, function(i) mean(x1[i,]))
	v1 <- sapply(1:samples, function(i) var(x1[i,]))

	list(x1 = x1,  m1=m1, v1=v1, parameters = pars, Xo = pars[1], time=times)
}




# Pars = {x, y, bx, by, dx, dy, cx, cy, K} 
crowley_ibm <- function(Xo = c(500,4500), 
	parameters=c(0.11, .6, .1, .1, .1, 4, 10000), 
	times = seq(0,100,length=50), reps=1){
	samples <- length(times)
	N <- reps*samples
	maxtime <- max(times)
	pars <- c(Xo, parameters)
	o <- .C("crowley", double(N), double(N), as.double(pars), as.integer(samples), as.integer(reps), as.double(maxtime) )

	x1 = matrix(o[[1]], samples, reps)
	m1 <- sapply(1:samples, function(i) mean(x1[i,]))
	v1 <- sapply(1:samples, function(i) var(x1[i,]))

	x2 = matrix(o[[2]], samples, reps)
	m2 <- sapply(1:samples, function(i) mean(x2[i,]))
	v2 <- sapply(1:samples, function(i) var(x2[i,]))

	list(x1 = x1, x2 = x2, m1=m1, m2=m2, v1=v1, v2=v2, parameters = parameters, Xo = Xo)
}


# Pars = {x, y, bx, by, dx, dy, cx, cy, K} 
metapop_ibm <- function(Xo = c(500,500), parameters=c(0.2, .6, .1, .1, .1, .1, 10000), times = seq(0,100,length=50), reps=1){
	samples <- length(times)
	N <- reps*samples
	maxtime <- max(times)
	pars <- c(Xo, parameters)
	o <- .C("metapop", double(N), double(N), as.double(pars), as.integer(samples), as.integer(reps), as.double(maxtime) )

	x1 = matrix(o[[1]], samples, reps)
	m1 <- sapply(1:samples, function(i) mean(x1[i,]))
	v1 <- sapply(1:samples, function(i) var(x1[i,]))

	x2 = matrix(o[[2]], samples, reps)
	m2 <- sapply(1:samples, function(i) mean(x2[i,]))
	v2 <- sapply(1:samples, function(i) var(x2[i,]))

	list(x1 = x1, x2 = x2, m1=m1, m2=m2, v1=v1, v2=v2, parameters = parameters, Xo = Xo)
}



# Pars = {b, ue, ul, up, ua, ae, al, ap, cle, cap, cae, Vol} 
gamma_beetles_ibm <- function(Xo = c(100,0,0,0), 
						parameters= c(b=5, ue=0, ul=0.001, up=0, ua=0.01, ae=1.3, al=.1, ap=1.5, cle=0.2, cap=0.1, cae=5, V=100),
						K = 10,
						times = seq(0,4000,length=400),
						reps = 1 ){
	samples <- length(times)
	N <- reps*samples
	maxtime <- max(times)
	## Note that the C code is using a very different convention for order of parameters
	pars <- c(parameters[c("ae", "al", "ap")], parameters[c("ue", "ul", "up", "ua")], parameters["b"], K, parameters[c("cle","cae","cap")], parameters["V"])
	n_rates = 6*K+2;
	n_states = 3*K+1;
	max_time = 200;

print(pars)

	# start at beginning of each age class
	inits = integer(n_states)
	inits[1] = Xo[1]
	inits[K+1] = Xo[2]
	inits[2*K+1] = Xo[3]
	inits[3*K+1] = Xo[4]

	o <- .C("gamma_beetles",  
				as.integer(inits),
				as.double(pars), 
				as.integer(n_rates),
				as.integer(n_states),
				as.double(maxtime),
				as.integer(samples), 
				as.integer(reps), 
				integer(N), # 8  
				integer(N), # 9 
				integer(N), #10
				integer(N) )#11 

	calc_moments <- function(j){
		x <- matrix(o[[7+j]], samples, reps)
		m <- sapply(1:samples, function(i) mean(x[i,]))
		v <- sapply(1:samples, function(i) var(x[i,]))
		list(m=m,v=v)
	}
	moments <- sapply(1:4, calc_moments)

	list(E = o[[8]], L = o[[9]], P = o[[10]], A = o[[11]], mv=moments, parameters = parameters, Xo = o[[1]], times = times)
}









# Pars = {E, L, P, A, b, ue, ul, up, ua, ae, al, ap, cle, cap, cae} 
beetles_ibm <- function(Xo = c(100,0,0,0), 
						parameters= c(5., 0, 0.001, 0, 0.003, 1/3.8, 1/(20.2-3.8), 1/(25.5-20.2), 0.01, 0.004, 0.01, 100),
						times = seq(0,1000,length=500),
						reps = 1 ){
	samples <- length(times)
	N <- reps*samples
	maxtime <- max(times)
	pars <- c(Xo, parameters)
	o <- .C("beetles", double(N), double(N), double(N), double(N),  as.double(pars), as.integer(samples), as.integer(reps), as.double(maxtime) )

	calc_moments <- function(j){
		x <- matrix(o[[j]], samples, reps)
		m <- sapply(1:samples, function(i) mean(x[i,]))
		v <- sapply(1:samples, function(i) var(x[i,]))
		list(m=m,v=v)
	}
	moments <- sapply(1:4, calc_moments)

	list(E = o[[1]], L = o[[2]], P = o[[3]], A = o[[4]], mv=moments, parameters = parameters, Xo = Xo)
}
