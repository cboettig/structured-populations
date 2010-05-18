
require(odesolve)


# A one-dimensional test case: b(n) = rn(K-n) 
oneD <- function(t, y, p){
	# \dot x = \alpha_1(x)
	yd1 <- p["r"]*y[1] *(1 - y[1]/p["K"])
	# \dot \sigma^2 = - \partial_x \alpha_1(x) \sigma^2 + \alpha_2(x)
	yd2 <- 2*p["r"]*(1-2*y[1]/p["K"] ) * y[2] + 
		p["r"]*y[1] *(1 + y[1]/p["K"])

	list(c(yd1, yd2))
}

logistic_parameters <- c(K=100, r=1, d = .1)
times <- seq(0,10,by=.1)
yo <- c(10, 0)
o <- lsoda(yo, times, oneD, logistic_parameters)


linnoise <- function(t,y,p, birth, death, Jacobian, Transition=0){
# General form for arbitrary dimensions
# Args:
#	t - time
#	y - statespace, a vector of the d states, followed by d varancies and (d^2-d)/2 covariances
#	p - vector of parameters.  
#	b - user-supplied birth rates, fn of (t,y,p)
#	d - user-supplied death rates, fn of (t,y,p)
#	J - user-supplied Jabobian of f, also fn of (t,y,p), returns a matrix
#	T - user-supplied two-step transition rates, assumed 0 if not provided
# Returns:
#	rate of change in y, (in all mean states, varainces, and covariances), as can be used by lsode.   
# Notes: 
#	Both y and and p can use named values if the f,g,J,T functions are written to use them, but will be slower
#	Typically provide wrapper fn to pass model; i.e. eqns <- function(t,y,p){linnoise(t,y,p,f,g,J,T)}

	n <- length(y)				# d states, d variances, (d^2-d)/2 covariances
	D <- -3/2 + sqrt(9/4+2*n)	# dimension of the system
	yd <- numeric(n)
	
	TwoStep_rates <- Transition(t,y,p)
	net_transitions <- sapply(1:D, function(i){ sum(TwoStep_rates[i,]-TwoStep_rates[,i]) } )
	
	# macroscopics
	yd[1:d] <- birth(t,y,p) - death(t,y,p)  + net_transitions

	# build the variance-covariance matrix
	M = matrix(0, D,D)
	diag(M) = y[(D+1):(2*D) ]
	k <- 2*D+1
	for (i in 1:(D-1)){
		for (j in (i+1):D){
			M[i,j] <- y[k]
			k <- k+1
		}
	}
	
	E <- Jacobian(t,y,p) %*% M
	dM <- E + t(E) + diag( b(t,y,p)+d(t,y,p) ) - TwoStep

	# unfold the matrix, consider cholskey decomposition!
	yd[(D+1):(2*D) ] <- diag(dM) 
	k <- 2*D+1
	for (i in 1:(D-1)){
		for (j in (i+1):D){
			dM[i,j] -> yd[k]
			k <- k+1
		}
	}
	list(yd)
}

# Matrix of two-step transitions
T_crowley <- function(t,y,p){
	matrix(
		c(	0, 0, 
			0, 0),
		2,2)
}

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
	yd2 <- - p["d2"] * y[2] - p["c2"] * y[1] * y[2]
	c(yd1, yd2)
}


J_crowley <- function(t,y,p){
	j <- c(		( p["b1"]*(p["K"] - 2*y[1] - y[2]) - p["d1"] + p["c1"]*y[2] ),		(-p["b1"] * y[1]  + p["c1"] * y[1]  ), 
			-(p["b2"] + p["c2"])*y[2],												( p["b2"] * (p["K"] - y[1] - 2*y[2]) - p["d2"] - p["c2"] * y[1] ) ) 
	t(matrix(j,2,2))
}

# lsoda needs a function of only t,y,p


pop<-10000
crowley_parameters <- c(b1=.11/pop, b2=.6/pop, d1=0.1, d2=.1, c1=0.1/pop, c2=4/pop, K=pop)
#crowley_parameters <- c(b1=.2/pop, b2=.6/pop, d1=0.1, d2=.1, c1=0.1/pop, c2=.2/pop, K=pop)
times <- seq(0,4000,length=1000)
yo <- c(6, 5, 0, 0, 0) # (xo, yo, sigma_xo sigma_yo, cov)


#o <- lsoda(yo, times, twoD, crowley_parameters)
eqns <- function(t,y,p){ linnoise(t,y,p, b_crowley, d_crowley, J_crowley, T_crowley) }
crowley_sim <- lsoda(yo, times, eqns, crowley_parameters)

#png("crowley_noise.png")
par(mfrow=c(2,1))
plot(crowley_sim[,1], crowley_sim[,2], col="darkblue", lwd=3, type='l', ylim=c(0,8000), xlab="time",ylab="mean", cex.lab=1.3, main="Modified Crowley Model")
lines(crowley_sim[,1], crowley_sim[,3], col="darkgreen", lwd = 3)
plot(crowley_sim[,1], sqrt(crowley_sim[,4]), col="darkblue", lwd=3, type='l', ylim=c(0,8000), xlab="time",ylab="stdev", cex.lab=1.3 )
lines(crowley_sim[,1], sqrt(crowley_sim[,5]), col="darkgreen", lwd = 3)
#dev.off()



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
