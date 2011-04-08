# @file: ind_based_models.R
# @author: Carl Boettiger, <cboettig@gmail.com>
# @section DESCRIPTION wrapper to the C code containing the gillespie simulation for individual-based models.  


#   Pars = {n, e, a, K, h, i, Da,Dt} 
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
    pars=c("Xo" = 570, "e" = .5, "a" = 160, "K" = 1000, "h" = 200, "i" = 0, "Da" = 1, "Dt" = 100, "p"=2),
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



