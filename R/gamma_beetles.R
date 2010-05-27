
## Beetle model macroscopic eqns
b_beetles <- function(t,y,p){
	o <- numeric(40);
	o[1] <- p["b"]*sum(y[31:40])
	o
}

## Beetle model macroscopic eqns
d_beetles <- function(t,y,p){
	o <- numeric(40)
	o[1:10] <-  (p["ue"] + p["cle"]*sum(y[11:20]) + p["cae"]*sum(y[31:40]) )*y[1:10]  
	o[11:20] <- p["ul"]*y[11:20]
	o[21:30] <- (p["up"]+p["cap"]*sum(y[31:40]))*y[21:30] 
	o[31:40] <- p["ua"]*y[31:40]
	o
}
## Jacobian of f
J_beetles <- function(t,y,p){
	o <- matrix(0, 40, 40)
	d <- numeric(40)

	#diagonal
	d[1:10] <- -p["ue"] - p["ae"] -p["cle"]*sum(y[11:20]) - p["cae"]*sum(y[31:40])
	d[11:20] <- -p["ul"] - p["al"]
	d[21:30] <- -p["up"] - p["ap"] -p["cap"]*sum(y[31:40])
	d[31:39] <- -p["ua"] - p["ap"]
	d[40] <- -p["ua"] 
	diag(o) <- d

	# lower diagonal,
	o[matrix(c(2:11,1:10),10)] <-  p["ae"]
	o[matrix(c(12:21,11:20),10)] <- p["al"]
	o[matrix(c(22:31,21:30),10)] <- p["ap"]
	o[matrix(c(32:39,31:40),9)] <- p["ap"]

	# cannibalism
	o[matrix(c(1:10, 2:11),10)] <- -p["cle"]*y[1:10]
	o[matrix(c(1:10, 4:13),10)] <- -p["cae"]*y[1:10]
	o[1, 4] <- o[1,4] + p["b"]*sum(y[31:40])
	o[matrix(c(21:30, 24:33),10)] <- -p["cap"]*y[21:30]
	o
}




## Beetle model macroscopic eqns
T_beetles <- function(t,y,p){
	o <- matrix(0,40,40)
	o[matrix(c(1:10,2:11), 10, 2)] <-  p["ae"]*y[1:10]
	o[matrix(c(11:20,12:21),10)] <- p["al"]*y[11:20]
	o[matrix(c(21:30,22:31),10)] <- p["ap"]*y[21:30]
	o[matrix(c(31:39,32:40), 9)] <- p["ap"]*y[31:39]
	o
}

### No need to have structure in adults!

beetles_example <- function(){

	volume <- 100
	beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0.00001, ua = 0.01, 
						ae = 1.3, al = .1, ap = 1.5,
						cle = .2, cap = .1, cae = 5, V=volume)
	times <- seq(0,900,length=100)
	Xo <- numeric(40)
	Xo[1] <- 100


	beetle_data <- linear_noise_approx(Xo, times, beetle_pars, b_beetles, d_beetles, J_beetles, T_beetles, Omega=volume)
#	beetle_data <- mean_only(Xo, times, beetle_pars, b_beetles, d_beetles, J_beetles, T_beetles, Omega=volume)

	data <- sapply(1:8, function(i){  rowSums(beetle_data[,(2+10*(i-1)):(1+10*i)] ) } )
	cols = c("yellow", "yellowgreen", "lightgreen", "darkgreen");

#	par(mfrow=c(2,1))


	m<- max(data[,1:4])
	plot(times, data[,1], col=cols[1], lwd=3, type='l', ylim=c(0,m), ylab="means", main="Gamma Aging" )
	for (i in 2:4) {
		lines(times, data[,i], col=cols[i], lwd=3)
	}
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )

	v<- max(sqrt(data[,5:6]))
	plot(times, sqrt(data[,5]), col=cols[1], lwd=3, type='l', ylim=c(0,v) )
	for (i in 2:2) {
		lines(times, sqrt(data[,4+i]), col=cols[i], lwd=3)
	}



}
