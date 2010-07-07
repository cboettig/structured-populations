## Simulate beetle dynamics


k <-10
eggs <- 1:k	
larva <- (k+1):(2*k)
pupa <- (2*k+1):(3*k)
adults <- 3*k+1

## Beetle model macroscopic eqns
b_gamma <- function(t,y,p){
	o <- numeric(adults);
	o[1] <- p["b"]*y[adults]
	o
}

## Beetle model macroscopic eqns
d_gamma <- function(t,y,p){
	o <- numeric(adults)

	o[eggs] <-  (p["ue"] + p["cle"]*sum(y[larva]) + p["cae"]*y[adults])*y[eggs]  
	o[larva] <- p["ul"]*y[larva]
	o[pupa] <- (p["up"]+p["cap"]*y[adults])*y[pupa] 
	o[adults] <- p["ua"]*y[adults]
	o
}
## Jacobian of f
J_gamma <- function(t,y,p){
	o <- matrix(0, adults, adults)
	d <- numeric(adults)

	#diagonal
	d[eggs] <- -p["ue"] - p["ae"] -p["cle"]*sum(y[larva]) - p["cae"]*y[adults]
	d[larva] <- -p["ul"] - p["al"]
	d[pupa] <- -p["up"] - p["ap"] -p["cap"]*y[adults]
	d[adults] <- -p["ua"] 
	diag(o) <- d

	# lower diagonal,
	o[matrix(c(eggs+1,eggs),k)] <-  p["ae"]
	o[matrix(c(larva+1,larva),k)] <- p["al"]
	o[matrix(c(pupa+1,pupa),k)] <- p["ap"]

	# cannibalism
	o[matrix(c(eggs, eggs+k),k)] <- -p["cle"]*y[eggs]
	o[matrix(c(eggs, rep(adults,k)),k)] <- -p["cae"]*y[eggs]
	o[1, adults] <- o[1,adults] + p["b"]*y[adults]
	o[matrix(c(pupa, rep(adults,k)),k)] <- -p["cap"]*y[pupa]
	o
}




## Beetle model macroscopic eqns
T_gamma <- function(t,y,p){
	o <- matrix(0,adults,adults)
	o[matrix(c(eggs,eggs+1),  k)] <-  p["ae"]*y[eggs]
	o[matrix(c(larva,larva+1),k)] <- p["al"]*y[larva]
	o[matrix(c(pupa,pupa+1), k)] <- p["ap"]*y[pupa]
	o
}

### No need to have structure in adults!

gamma_example <- function(){



	volume <- 100
#	beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0, ua = 0.01, 
#						ae = .13*k, al = .01*k, ap = .15*k,
#						cle = .2, cap = .1, cae = 5, V=volume)
## 
	beetle_pars <- c(	b=5, ue= 0, ul = 0.001, up = 0, ua = .001, 
						ae = .1*k, al = .01*k, ap = .1*k,
						cle = 1, cap = .4, cae = 1, V=volume)

	times <- seq(0,400,length=100)
	Xo <- numeric(adults)
	Xo[1] <- 250
	Xo[11] <- 50
	Xo[21] <- 10
	Xo[31] <- 50
#	Xo = rep(10, adults)

	beetle_data <- linear_noise_approx	(Xo, times, beetle_pars, b_gamma, d_gamma, J_gamma, T_gamma, Omega=volume)
	ibm <- gamma_beetles_ibm(Xo=c(Xo[1], Xo[11], Xo[21], Xo[31]), par=beetle_pars, time = times, reps = 100)


	save(list=ls(), file = "gamma_beetles_data.Rdat")

	# Collapse the pseudo-classes
	data <- matrix(0, length(times), 8)
	for( i in 1:3) {
		data[,i] = rowSums(beetle_data[,(2+k*(i-1)):(1+k*i)] )
	}
	data[,4] = beetle_data[,3*k+2]
	vars_start <- 3*k+3
	for( i in 1:3) {
		data[,4+i] = rowSums(beetle_data[,(vars_start+k*(i-1)):(vars_start+1+k*i)] )
	}
	data[,8] = beetle_data[,1+2*(3*k+1)]


	png("poisson_noise.png", 800, 800)
	# Plot results
	cols = c("yellow", "yellowgreen", "lightgreen", "darkgreen");
	par(mfrow=c(2,1))
	m<- max(data[,1:4])
	plot(times, data[,1], col=cols[1], lwd=3, type='l', ylim=c(0,m), ylab="means", main="Gamma Aging" )
	for (i in 2:4) {
		lines(times, data[,i], col=cols[i], lwd=3)
	}
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )
	
	points(times, ibm$mv[[1,1]], col="yellow")	
	points(times, ibm$mv[[1,2]], col="yellowgreen")	
	points(times, ibm$mv[[1,3]], col="lightgreen")	
	points(times, ibm$mv[[1,4]], col="darkgreen")	



	v<- max(data[,5:8])
	plot(times, (data[,5]), col=cols[1], lwd=3, type='l', ylim=c(0,v) )
	for (i in 2:4) {
		lines(times, (data[,4+i]), col=cols[i], lwd=3)
	}


	points(beetle_data[,1], sqrt(ibm$mv[[2,1]]), col="yellow")	
	points(beetle_data[,1], sqrt(ibm$mv[[2,2]]), col="yellowgreen")	
	points(beetle_data[,1], sqrt(ibm$mv[[2,3]]), col="lightgreen")	
	points(beetle_data[,1], sqrt(ibm$mv[[2,4]]), col="darkgreen")

	dev.off()

}
