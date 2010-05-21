# @file: ind_based_models.R
# @author: Carl Boettiger, <cboettig@gmail.com>
# @section DESCRIPTION wrapper to the C code containing the gillespie simulation for individual-based models.  

# Pars = {x, y, bx, by, dx, dy, cx, cy, K} 
crowley_ibm <- function(Xo = c(500,4500), parameters=c(0.11/10000, .6/10000, .1, .1, .1/10000, 4/10000, 10000), times = seq(0,1000,length=500), reps=1){
	samples <- length(times)
	N <- reps*(1+samples)
	maxtime <- max(times)
	pars <- c(Xo, parameters)
	o <- .C("crowley", double(N), double(N), as.double(pars), as.integer(samples), as.integer(reps), as.double(maxtime) )
	list(x1 = o[[1]], x2 = o[[2]], parameters = parameters, Xo = Xo)
}




# Pars = {E, L, P, A, b, ue, ul, up, ua, ae, al, ap, cle, cap, cae} 
beetles_ibm <- function(Xo = c(100,0,0,0), 
						parameters= c(5., 0, 0.001, 0, 0.003, 1/3.8, 1/(20.2-3.8), 1/(25.5-20.2), 0.01, 0.004, 0.01),
						times = seq(0,1000,length=500)
						){
	N <- length(times)
	maxtime <- max(times)
	pars <- c(Xo, parameters)
	o <- .C("beetles", double(N), double(N), double(N), double(N),  as.double(pars), as.integer(N), as.double(maxtime) )
	list(E = o[[1]], L = o[[2]], P = o[[3]], A = o[[4]], parameters = parameters, Xo = Xo)
}
