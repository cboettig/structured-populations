
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
						ae = .01, al = .01, ap = .1,
						cle = 3, cap = .4, cae = 1, V=volume)
	times <- seq(0,4000,length=100)
	Xo <- c(100,0,0,0)
	beetle_data <- linear_noise_approx(Xo, times, beetle_pars, b_beetles, d_beetles, J_beetles, T_beetles, Omega=volume) 
#	ibm <- beetles_ibm(Xo=Xo, par=beetle_pars, time=times, reps=40)


	

#	png("beetles_bugs.png")
	par(mfrow=c(2,2))
	m <- max(beetle_data[,2:5])*1.2
	plot(beetle_data[,1], beetle_data[,2], type = 'l', col="yellow", 
		lwd=3, ylim=c(0,m), xlab="time", ylab="mean", cex.lab=1.3, main="Beetle ELPA model" )
	lines(beetle_data[,1], beetle_data[,3], col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], beetle_data[,4], col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], beetle_data[,5], col="darkgreen", lwd = 3)	
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )

#	points(beetle_data[,1], ibm$mv[[1,1]], col="yellow")	
#	points(beetle_data[,1], ibm$mv[[1,2]], col="yellowgreen")	
#	points(beetle_data[,1], ibm$mv[[1,3]], col="lightgreen")	
#	points(beetle_data[,1], ibm$mv[[1,4]], col="darkgreen")	


	v <- max(sqrt(beetle_data[,6:9]))
	plot(beetle_data[,1], sqrt(beetle_data[,6]), type = 'l', col="yellow",
		lwd=3, ylim=c(0,v), xlab="time", ylab="stdev", cex=1.3 )
	lines(beetle_data[,1], sqrt(beetle_data[,7]), col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], sqrt(beetle_data[,8]), col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], sqrt(beetle_data[,9]), col="darkgreen", lwd = 3)

	eps <- .1
	v <- max(sqrt(beetle_data[,6:9])/(eps+beetle_data[,2:5]))
	plot(beetle_data[,1], sqrt(beetle_data[,6])/beetle_data[,2], type = 'l', col="yellow",
		lwd=3, ylim=c(0,v), xlab="time", ylab="stdev/mean", cex=1.3 )
	lines(beetle_data[,1], sqrt(beetle_data[,7])/beetle_data[,3], col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], sqrt(beetle_data[,8])/beetle_data[,4], col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], sqrt(beetle_data[,9])/beetle_data[,5], col="darkgreen", lwd = 3)


#	points(beetle_data[,1], sqrt(ibm$mv[[2,1]]), col="yellow")	
#	points(beetle_data[,1], sqrt(ibm$mv[[2,2]]), col="yellowgreen")	
#	points(beetle_data[,1], sqrt(ibm$mv[[2,3]]), col="lightgreen")	
#	points(beetle_data[,1], sqrt(ibm$mv[[2,4]]), col="darkgreen")

#	dev.off()


	vu <- max(beetle_data[,10:15])
	vl <- min(beetle_data[,10:15])
	plot(beetle_data[,1], beetle_data[,10], type = 'l', col="yellow",
		lwd=3, ylim=c(vl,vu), xlab="time", ylab="cov", cex=1.3 )
	lines(beetle_data[,1], beetle_data[,11], col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], beetle_data[,12], col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], beetle_data[,13], col="darkgreen", lwd = 3)
	lines(beetle_data[,1], beetle_data[,14], col="lightblue", lwd = 3)
	lines(beetle_data[,1], beetle_data[,15], col="blue", lwd = 3)




#	dev.off()
}
