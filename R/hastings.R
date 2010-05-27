b_h <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"]*y[1]*(1 - y[1] - y[2])  
	# \dot y = \beta_1(x,y)
	yd2 <-0 
	c(yd1, yd2)
}

d_h <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["d1"] * y[1] 
	# \dot y = \beta_1(x,y)
	yd2 <- p["d2"] * y[2] 
	c(yd1, yd2)
}


J_h <- function(t,y,p){
	t(matrix( c(		( p["b1"]*(1 - 2*y[1] - y[2]) - p["d1"]  ),
				-p["b1"]*y[1], 
				0,
				0 
			),2,2))
}

# Matrix of two-step transitions
T_h <- function(t,y,p){
	t(matrix(
		c(	0, p["c1"]*y[1]*y[2], 0, 0),
		2,2))
}


hastings_example <- function(){

	pop <- 10000
	parameters <- c(b1=.2, b2=0, 
							d1=0.1, d2=.1, c1=0.5, 
							c2=0, K=pop)

	times <- seq(0,200,length=50)
	Xo <- c(1500, 1500)
	sim <- linear_noise_approx(
				Xo, times, parameters, 
				b_h, d_h, J_h,
				T_h, Omega=pop)
#	ibm <- metapop_ibm(Xo = Xo, parameters=crowley_parameters, times=times,reps=80)
	meta_parameters <- c(b1=.0, b2=.2, 
							d1=0.1, d2=.1, c1=0.5, 
							c2=0.0, K=pop)

	m_sim <- linear_noise_approx(
				Xo, times, meta_parameters, 
				b_meta, d_meta, J_meta,
				T_meta, Omega=pop)

	par(mfrow=c(2,1))
	m <- max(sim[,2:3])*1.4
	plot(sim[,1], sim[,2], col="darkblue", lwd=3, type='l', xlab="time",ylab="mean", ylim = c(0,m), cex.lab=1.3, main="Two-step vs one-step Noise")
	lines(sim[,1],sim[,3], col="darkgreen", lwd = 3)

	v <- max(sqrt(sim[,4:5]))*1.2
	plot(sim[,1], sqrt(sim[,4]), ylim = c(0,v), col="darkblue", lwd=3, type='l', xlab="time",ylab="stdev", cex.lab=1.3 )
	lines(sim[,1], sqrt(sim[,5]), col="darkgreen", lwd = 3)



}

