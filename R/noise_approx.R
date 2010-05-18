
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

# Transition=function(){0}
linnoise <- function(t,y,p, birth, death, Jacobian, Transition){
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
	
	transition_matrix <- Transition(t,y,p)
	net_transitions <- sapply(1:D, function(i){ sum(transition_matrix[,i]-transition_matrix[i,]) } )
	
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
	
	#build T matrix
	T <- transition_matrix + t(transition_matrix)
	diag(T) = -rowSums(T)

	E <- Jacobian(t,y,p) %*% M
	dM <- E + t(E) + diag( birth(t,y,p)+death(t,y,p) ) - T

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




