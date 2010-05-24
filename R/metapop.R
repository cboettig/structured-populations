b_meta <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"]*y[1]*(1 - y[1] - y[2])  
	# \dot y = \beta_1(x,y)
	yd2 <- p["b2"]*y[2]*(1 - y[1] - y[2]) 
	c(yd1, yd2)
}

d_meta <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["d1"] * y[1] 
	# \dot y = \beta_1(x,y)
	yd2 <- p["d2"] * y[2] 
	c(yd1, yd2)
}


J_meta <- function(t,y,p){
	t(matrix( c(		( p["b1"]*(1 - 2*y[1] - y[2]) - p["d1"]  ),
				-p["b1"]*y[1], 
				-p["b2"]*y[2],
				( p["b2"] * (1 - y[1] - 2*y[2]) - p["d2"] ) 
			),2,2))
}

# Matrix of two-step transitions
T_meta <- function(t,y,p){
	t(matrix(
		c(	0, 0, 
		p["c1"]*y[1]*y[2], 0),
		2,2))
}


metapop_example <- function(){

	pop <- 10000
	crowley_parameters <- c(b1=.2, b2=.6, 
							d1=0.1, d2=.1, c1=0.15, 
							c2=0.15, K=pop)

	times <- seq(0,100,length=50)
	Xo <- c(1500, 1500)
	meta_sim <- linear_noise_approx(
				Xo, times, crowley_parameters, 
				b_meta, d_meta, J_meta,
				T_meta, Omega=pop)

#	source("crowley.R")
	crowley_sim <- linear_noise_approx(
				Xo, times, crowley_parameters, 
				b_crowley, d_crowley, J_crowley,
				T_crowley, Omega=pop)

	ibm <- metapop_ibm(Xo = Xo, parameters=crowley_parameters, times=times,reps=80)
	ibm2 <- crowley_ibm(Xo = Xo, parameters=crowley_parameters, times=times,reps=80)


	png("crowley_noise2.png", 800,800)
	par(mfrow=c(2,1))
	m <- max(crowley_sim[,2:3])*1.4
	plot(crowley_sim[,1], crowley_sim[,2], col="darkblue", lwd=3, type='l', xlab="time",ylab="mean", ylim = c(0,m), cex.lab=1.3, main="Two-step vs one-step Noise")
	lines(crowley_sim[,1], crowley_sim[,3], col="darkgreen", lwd = 3)
	lines(crowley_sim[,1], meta_sim[,2], col="lightblue", lwd = 3, lty="dashed")
	lines(crowley_sim[,1], meta_sim[,3], col="lightgreen", lwd = 3, lty="dashed")
	legend("right", c("crolwey sim", "metapop sim"), pch=c('+', 'o'))


	points(crowley_sim[,1], ibm$m1, col="darkblue", pch='+')
	points(crowley_sim[,1], ibm$m2, col="darkgreen", pch='+')
	points(crowley_sim[,1], ibm2$m1, col="darkblue" )
	points(crowley_sim[,1], ibm2$m2, col="darkgreen")

		

	v <- max(sqrt(crowley_sim[,4:5]))*1.2
	plot(crowley_sim[,1], sqrt(crowley_sim[,4]), ylim = c(0,v), col="darkblue", lwd=3, type='l', xlab="time",ylab="stdev", cex.lab=1.3 )
	lines(crowley_sim[,1], sqrt(crowley_sim[,5]), col="darkgreen", lwd = 3)
	lines(crowley_sim[,1], sqrt(meta_sim[,4]), col="lightblue", lwd = 4, lty="dotted")
	lines(crowley_sim[,1], sqrt(meta_sim[,5]), col="lightgreen", lwd = 4, lty="dotted")

	points(crowley_sim[,1], sqrt(ibm$v1), col="darkblue", pch='+')
	points(crowley_sim[,1], sqrt(ibm$v2), col="darkgreen", pch='+')
	points(crowley_sim[,1], sqrt(ibm2$v1), col="darkblue" )
	points(crowley_sim[,1], sqrt(ibm2$v2), col="darkgreen")
	legend("bottomright", c("crowley x", "crowley y", "metapop x", "metapop y"), lty=c(1,1,2,2), col=c("darkblue", "darkgreen", "lightblue", "lightgreen") )

	dev.off()
}

