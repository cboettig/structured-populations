require(stochPop)
source("../R/onestep_beetles.R")
source("../R/plots.R")

volume <- 100
beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0, ua = .001, 
						ae = .1, al = .01, ap = .1,
						cle = 1, cap = .4, cae = 1, V=volume)
times <- seq(0,500,length=100)
Xo <- c(100,0,0,0)

ibm <- beetles_ibm(Xo=Xo, par=beetle_pars, time=times, reps=100)

beetle_data_one <- linear_noise_approx(Xo, times, beetle_pars, 
					b_beetles_one, d_beetles_one, J_beetles_one,
					T_beetles_one, Omega=volume) 
beetle_data <- linear_noise_approx(Xo, times, beetle_pars,
					b_beetles, d_beetles, J_beetles,
					T_beetles, Omega=volume) 

plot_stdev(beetle_data_one)
plot_stdev(beetle_data, overlay=T, ibm=ibm)




source("../R/gamma_beetles.R")
k <- 10; adults <- 3*k+1
gamma_beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0, ua = .001, 
						ae = .1*k, al = .01*k, ap = .1*k,
						cle = 1, cap = .4, cae = 1, V=volume)
Xo <- numeric(adults)
Xo[1] <- 100
beetle_data_gamma <- linear_noise_approx(Xo, times, gamma_beetle_pars,
					b_gamma, d_gamma, J_gamma, T_gamma, Omega=volume)
ibm_gamma <- gamma_beetles_ibm(Xo=c(Xo[1], Xo[11], Xo[21], Xo[31]),
					par=beetle_pars, time = times, reps = 100)
beetle_gamma <- collapse_gamma_classes(beetle_data_gamma, k)


plot_stdev(beetle_gamma, ibm=ibm_gamma)





