
## Beetle model macroscopic eqns
b_beetles <- function(t,y,p){
	f1 <- p["b"]*y[4] 
	f2 <- 0 
	f3 <- 0
	f4 <- 0
	c(f1, f2, f3, f4)
}

## Beetle model macroscopic eqns
d_beetles <- function(t,y,p){
	f1 <-  (p["ue"] + p["cle"]*y[2] + p["cae"]*y[4] )*y[1]  
	f2 <- p["ul"]*y[2]
	f3 <- (p["up"]+p["cap"]*y[4])*y[3] 
	f4 <- p["ua"]*y[4] 
	c(f1, f2, f3, f4)
}

## Jacobian of f
J_beetles <- function(t,y,p){
 j <- c(	-(p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4]),
			-p["cle"]*y[1],
			0, 
			p["b"] - p["cae"]*y[1],

			p["ae"],
			-p["ul"]-p["al"], 
			0, 
			0,

			0, 
			p["al"], 
			-p["up"]-p["cap"]*y[4]-p["ap"], 
			-p["cap"]*y[3],

			0,
			0,
			p["ap"],
			-p["ua"] )
 t(matrix(j,4,4))
}




## Beetle model macroscopic eqns
T_beetles <- function(t,y,p){
	t(matrix(c(
		0,		p["ae"]*y[1],	0,				0,
		0,		0,				p["al"]*y[2],	0,
		0,		0,				0,				p["ap"]*y[3],
		0,		0,				0,				0
		), 4,4))
}



beetles_example <- function(){

	volume <- 100
	beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0, ua = .001, 
						ae = .1, al = .01, ap = .1,
						cle = 1, cap = .4, cae = 1, V=volume)
	times <- seq(0,500,length=100)
	Xo <- c(100,0,0,0)
	beetle_data <- linear_noise_approx(Xo, times, beetle_pars, b_beetles, d_beetles, J_beetles, T_beetles, Omega=volume) 
#	ibm <- beetles_ibm(Xo=Xo, par=beetle_pars, time=times, reps=40)

#	png("beetles_bugs.png")
	par(mfrow=c(2,1))
	plot_means(beetle_data)
	plot_stdev(beetle_data)
	plot_covar(beetle_data)
#	dev.off()
}
