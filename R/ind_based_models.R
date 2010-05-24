# @file: ind_based_models.R
# @author: Carl Boettiger, <cboettig@gmail.com>
# @section DESCRIPTION wrapper to the C code containing the gillespie simulation for individual-based models.  

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
