
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
 t( matrix(c(
	-(p["ue"]+p["cle"]*y[2]+p["cae"]*y[4]),			-p["cle"]*y[1],	0,	p["b"] - p["cae"]*y[1],
	  0,			-p["ul"],	0,						0,
	  0,				0,	-p["up"]-p["cap"]*y[4],		-p["cap"]*y[3],
	  0,				0,		0,						-p["ua"] 
	  ), 4, 4))

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
	beetle_pars <- c(	b=5, ue= 0, ul = .001, up = 0, ua = .003, 
						ae = 1/3.8, al = 1/(20.2-3.8), ap = 1/(25.5-20.2),
						cle = 0.1, cap = 0.04, cae = 0.1, V=volume)
	times <- seq(0,200,length=50)
	Xo <- c(1000,6000,10,2000)
	beetle_data <- linear_noise_approx(Xo, times, beetle_pars, b_beetles, d_beetles, J_beetles, T_beetles, Omega=volume) 
	ibm <- beetles_ibm(Xo=Xo, par=beetle_pars, time=times, reps=4)


	

	png("beetles_bugs.png")
	par(mfrow=c(2,1))
	m <- max(beetle_data[,2:5])*1.2
	plot(beetle_data[,1], beetle_data[,2], type = 'l', col="yellow", 
		lwd=3, ylim=c(0,m), xlab="time", ylab="mean", cex.lab=1.3, main="Beetle ELPA model" )
	lines(beetle_data[,1], beetle_data[,3], col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], beetle_data[,4], col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], beetle_data[,5], col="darkgreen", lwd = 3)	
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )

	points(beetle_data[,1], ibm$mv[[1,1]], col="yellow")	
	points(beetle_data[,1], ibm$mv[[1,2]], col="yellowgreen")	
	points(beetle_data[,1], ibm$mv[[1,3]], col="lightgreen")	
	points(beetle_data[,1], ibm$mv[[1,4]], col="darkgreen")	


	v <- max(sqrt(beetle_data[,6:9]))
	plot(beetle_data[,1], sqrt(beetle_data[,6]), type = 'l', col="yellow",
		lwd=3, ylim=c(0,v), xlab="time", ylab="stdev", cex=1.3 )
	lines(beetle_data[,1], sqrt(beetle_data[,7]), col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], sqrt(beetle_data[,8]), col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], sqrt(beetle_data[,9]), col="darkgreen", lwd = 3)

	points(beetle_data[,1], sqrt(ibm$mv[[2,1]]), col="yellow")	
	points(beetle_data[,1], sqrt(ibm$mv[[2,2]]), col="yellowgreen")	
	points(beetle_data[,1], sqrt(ibm$mv[[2,3]]), col="lightgreen")	
	points(beetle_data[,1], sqrt(ibm$mv[[2,4]]), col="darkgreen")	
	dev.off()
}
