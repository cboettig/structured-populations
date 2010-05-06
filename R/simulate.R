# simulate.R
beetle_sim <- function(	state = c(100, 10, 10, 10, 10), 
						pars = c(	5,
									0,
									0.001,
									0.0,
									0.003,
									3.8,
									3.8+16.4, 
									3.8+16.4+5, 
									0.01, 
									0.004, 
									0.01, 
									3.8+8
								), 
						dt = 14,
						T = 200,
						seed = 123 ){
	out <- .C("_Z10beetle_simPiPdS0_S0_S_", as.integer(state), as.double(pars), as.double(dt), as.double(T), as.integer(seed) )
	list(state=out[[1]], pars = out[[2]], dt = dt)
}


ensemble <- function(	state = c(450, 0, 0, 0, 30), 
						initial = state,
						pars = c(	5,
									0,
									0.001,
									0.0,
									0.003,
									3.8,
									3.8+16.4, 
									3.8+16.4+5, 
									0.01, 
									0.004, 
									0.01, 
									3.8+8
						), 
						dt = 14, 
						seed = 123,
						reps = 10, 
						probs = double(5),
						nstates = 5)
{
	out <- .C(	"_Z8ensemblePiS_PdS0_S_S_S0_S_", 
				as.integer(state), 
				as.integer(initial), 
				as.double(pars), 
				as.double(dt), 
				as.integer(seed), 
				as.integer(reps), 
				as.double(probs), 
				as.integer(nstates) )
	list(probs=out[[7]], pars = out[[2]], dt = dt)
}



