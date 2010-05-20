b_crowley <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"]*y[1]*(p["K"] - y[1] - y[2])  + p["c1"]*y[1]*y[2]
	# \dot y = \beta_1(x,y)
	yd2 <- p["b2"]*y[2]*(p["K"] - y[1] - y[2]) 
	c(yd1, yd2)
}

d_crowley <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["d1"] * y[1] 
	# \dot y = \beta_1(x,y)
	yd2 <- p["d2"] * y[2] + p["c2"] * y[1] * y[2]
	c(yd1, yd2)
}


J_crowley <- function(t,y,p){
	j <- c(		( p["b1"]*(p["K"] - 2*y[1] - y[2]) - p["d1"] + p["c1"]*y[2] ),		(-p["b1"] * y[1]  + p["c1"] * y[1]  ), 
			-(p["b2"] + p["c2"])*y[2],												( p["b2"] * (p["K"] - y[1] - 2*y[2]) - p["d2"] - p["c2"] * y[1] ) ) 
	t(matrix(j,2,2))
}

# Matrix of two-step transitions
T_crowley <- function(t,y,p){
	matrix(
		c(	0, 0, 
			0, 0),
		2,2)
}


crowley_example <- function(){
	pop <- 10000
	crowley_parameters <- c(b1=.11/pop, b2=.6/pop, d1=0.1, d2=.1, c1=0.1/pop, c2=4/pop, K=pop)
	#crowley_parameters <- c(b1=.2/pop, b2=.6/pop, d1=0.1, d2=.1, c1=0.1/pop, c2=.2/pop, K=pop)
	times <- seq(0,2000,length=1000)
	Xo <- c(595, 4550)

	crowley_sim <-linear_noise_approx(Xo, times, crowley_parameters, 
										b_crowley, d_crowley, J_crowley,
										T_crowley, Omega=pop)

	ibm <- crowley_ibm(Xo = Xo, par=crowley_parameters, time=times)

	#png("crowley_noise.png")
	par(mfrow=c(2,1))
	m <- max(crowley_sim[,2:3])*1.2
	plot(crowley_sim[,1], crowley_sim[,2], col="darkblue", lwd=3, type='l', xlab="time",ylab="mean", ylim = c(0,m), cex.lab=1.3, main="Modified Crowley Model")
	lines(crowley_sim[,1], crowley_sim[,3], col="darkgreen", lwd = 3)
	legend("right", c("competitor", "colonist"), lty=1, col=c("darkblue", "darkgreen") )

		
	points(crowley_sim[,1], ibm$x1, col="darkblue")
	points(crowley_sim[,1], ibm$x2, col="darkgreen")

	print(c(sd(ibm$x1), sd(ibm$x2)))	

	v <- max(sqrt(crowley_sim[,4:5]))
	plot(crowley_sim[,1], sqrt(crowley_sim[,4]), ylim = c(0,v), col="darkblue", lwd=3, type='l', xlab="time",ylab="stdev", cex.lab=1.3 )
	lines(crowley_sim[,1], sqrt(crowley_sim[,5]), col="darkgreen", lwd = 3)
	#dev.off()
}

