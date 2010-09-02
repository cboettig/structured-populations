require(stochPop)
source("../R/onestep_beetles.R")

	volume <- 100
	beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0, ua = .001, 
						ae = .1, al = .01, ap = .1,
						cle = 1, cap = .4, cae = 1, V=volume)
	times <- seq(0,500,length=100)
	Xo <- c(100,0,0,0)
	beetle_data_one <- linear_noise_approx(Xo, times, beetle_pars, b_beetles_one, d_beetles_one, J_beetles_one, T_beetles_one, Omega=volume) 
	beetle_data <- linear_noise_approx(Xo, times, beetle_pars, b_beetles, d_beetles, J_beetles, T_beetles, Omega=volume) 

#	ibm <- beetles_ibm(Xo=Xo, par=beetle_pars, time=times, reps=40)

	plot_means(beetle_data)


	plot_stdev(beetle_data)
	plot_covar(beetle_data)







