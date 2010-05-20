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
	pop<-10000
	crowley_parameters <- c(b1=.11/pop, b2=.6/pop, d1=0.1, d2=.1, c1=0.1/pop, c2=4/pop, K=pop)
	#crowley_parameters <- c(b1=.2/pop, b2=.6/pop, d1=0.1, d2=.1, c1=0.1/pop, c2=.2/pop, K=pop)
	times <- seq(0,4000,length=1000)
	yo <- c(6, 5, 0, 0, 0) # (xo, yo, sigma_xo sigma_yo, cov)
	eqns <- function(t,y,p){ linnoise(t,y,p, b_crowley, d_crowley, J_crowley, T_crowley) }
	crowley_sim <- lsoda(yo, times, eqns, crowley_parameters)
	#png("crowley_noise.png")
	par(mfrow=c(2,1))
	plot(crowley_sim[,1], crowley_sim[,3], col="darkblue", lwd=3, type='l', xlab="time",ylab="mean", cex.lab=1.3, main="Modified Crowley Model")
	lines(crowley_sim[,1], crowley_sim[,2], col="darkgreen", lwd = 3)
	plot(crowley_sim[,1], sqrt(crowley_sim[,5]), col="darkblue", lwd=3, type='l', xlab="time",ylab="stdev", cex.lab=1.3 )
	lines(crowley_sim[,1], sqrt(crowley_sim[,4]), col="darkgreen", lwd = 3)
	#dev.off()
}

