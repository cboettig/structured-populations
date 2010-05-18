# y is given as  c( state, variances, covariances )  states = d, then length = d+d+(d^2-d)/2 = 1.5*d+d^2/2 = n
#d <- 4
#n <- 1.5*d+d^2/2



## Beetle model parameters
beetle_pars <- c(	b=5, ue= 0, ul = .001, up = 0, ua = .03, ae = 1/3.8, al = 1/(20.2-3.8), ap = 1/(25.5-20.2), cle = 0.01, cap = 0.004, cae = 0.01)

## Beetle model macroscopic eqns
f_beetles <- function(t,y,p){
	f1 <- p["b"]*y[4]  - (p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4] )*y[1]  
	f2 <- p["ae"]*y[1] - (p["ul"] + p["al"])*y[2]
	f3 <- p["al"]*y[2] - (p["up"] + p["ap"])*y[3] 
	f4 <- p["ap"]*y[3] - p["ua"]*y[4] 
	c(f1, f2, f3, f4)
}


## Beetles second jump moment
g_beetles <- function(t,y,p){
	g1 <- p["b"]*y[4]  + (p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4] )*y[1]  
	g2 <- p["ae"]*y[1] + (p["ul"] + p["al"])*y[2]
	g3 <- p["al"]*y[2] + (p["up"] + p["ap"])*y[3] 
	g4 <- p["ap"]*y[3] + p["ua"]*y[4] 
	c(g1, g2, g3, g4)
}

## Jacobian of f
J_beetles <- function(t,y,p){
 j <- c( -(p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4]),	-p["cle"]*y[1],	0,			p["b"] - p["cae"]*y[1],
	  p["ae"],			-p["ul"]-p["al"],	0,			0,
	  0,				p["al"],	-p["up"]-p["ap"],	0,
	  0,				0,			p["ap"],	-p["ua"] )
 t(matrix(j,4,4))
}

## two-step transtion
T_beetles <- function(t,y,p){

}

times <- seq(0,1000,length=100)
d <- 4
n <- 1.5*d+d^2/2
yo_beetles <- numeric(n)
yo_beetles[1] = 100

beetle_eqns <- function(t,y,p){ linnoise(t,y,p, f_beetles, g_beetles, J_beetles) }
beetle_data <- lsoda(yo_beetles, times, beetle_eqns, beetle_pars)


png("beetles_noise.png")
par(mfrow=c(2,1))
m <- max(beetle_data[,2:5])
plot(beetle_data[,1], beetle_data[,2], type = 'l', col="yellow", lwd=3, ylim=c(0,m), xlab="time", ylab="mean", cex.lab=1.3, main="Beetle ELPA model" )
lines(beetle_data[,1], beetle_data[,3], col="yellowgreen", lwd=3)	
lines(beetle_data[,1], beetle_data[,4], col="lightgreen", lwd = 3)	
lines(beetle_data[,1], beetle_data[,5], col="darkgreen", lwd = 3)	

v <- max(sqrt(beetle_data[,6:9]))
plot(beetle_data[,1], sqrt(beetle_data[,6]), type = 'l', col="yellow", lwd=3, ylim=c(0,v), xlab="time", ylab="stdev", cex=1.3 )
lines(beetle_data[,1], sqrt(beetle_data[,7]), col="yellowgreen", lwd=3)	
lines(beetle_data[,1], sqrt(beetle_data[,8]), col="lightgreen", lwd = 3)	
lines(beetle_data[,1], sqrt(beetle_data[,9]), col="darkgreen", lwd = 3)	
dev.off()
