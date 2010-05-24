b_crowley <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"]*y[1]*(1 - y[1] - y[2])  + p["c1"]*y[1]*y[2]
	# \dot y = \beta_1(x,y)
	yd2 <- p["b2"]*y[2]*(1 - y[1] - y[2]) 
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
	t(matrix( c(		( p["b1"]*(1 - 2*y[1] - y[2]) - p["d1"] + p["c1"]*y[2] ),
				-p["b1"]*y[1]  + p["c1"]*y[1], 
				-p["b2"]*y[2] - p["c2"]*y[2],
				( p["b2"] * (1 - y[1] - 2*y[2]) - p["d2"] - p["c2"]*y[1] ) 
			),2,2))
}

# Matrix of two-step transitions
T_crowley <- function(t,y,p){
	matrix(
		c(	0, 0, 
			0, 0),
		2,2)
}


crowley_example <- function(){

	pop <- 1000
	crowley_parameters <- c(b1=.2, b2=.6, 
							d1=0.1, d2=.1, c1=0.1, 
							c2=.1, K=pop)

	times <- seq(0,200,length=50)
	Xo <- c(500, 4500)
	ibm <- crowley_ibm(Xo = Xo, parameters=crowley_parameters, times=times,reps=20)

	crowley_sim <- linear_noise_approx(
				Xo, times, crowley_parameters, 
				b_crowley, d_crowley, J_crowley,
				T_crowley, Omega=pop)

	#png("crowley_noise.png")
	par(mfrow=c(2,1))
	m <- max(crowley_sim[,2:3])*1.2
	plot(crowley_sim[,1], crowley_sim[,2], col="darkblue", lwd=3, type='l', xlab="time",ylab="mean", ylim = c(0,m), cex.lab=1.3, main="Modified Crowley Model")
	lines(crowley_sim[,1], crowley_sim[,3], col="darkgreen", lwd = 3)
	legend("right", c("competitor", "colonist"), lty=1, col=c("darkblue", "darkgreen") )

		
	points(crowley_sim[,1], ibm$m1, col="darkblue", pch='+')
	points(crowley_sim[,1], ibm$m2, col="darkgreen", pch='+')

	v <- max(sqrt(crowley_sim[,4:5]))*2
	plot(crowley_sim[,1], sqrt(crowley_sim[,4]), ylim = c(0,v), col="darkblue", lwd=3, type='l', xlab="time",ylab="stdev", cex.lab=1.3 )
	lines(crowley_sim[,1], sqrt(crowley_sim[,5]), col="darkgreen", lwd = 3)

	points(crowley_sim[,1], sqrt(ibm$v1), col="darkblue", pch='+')
	points(crowley_sim[,1], sqrt(ibm$v2), col="darkgreen", pch='+')

#	dev.off()
}

