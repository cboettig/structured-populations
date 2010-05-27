
require(odesolve)

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
	net_transitions <- sapply(1:D, 
		function(i){
			sum(transition_matrix[,i]-transition_matrix[i,]) 
		})
	
	# macroscopics
	yd[1:D] <- birth(t,y,p) - death(t,y,p)  + net_transitions

	# build the variance-covariance matrix
	M = matrix(0, D,D)
	diag(M) = y[(D+1):(2*D) ]
	k <- 2*D+1
	for (i in 1:(D-1)){
		for (j in (i+1):D){
			M[i,j] <- y[k]
			M[j,i] <- y[k]
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

linear_noise_approx <- function(Xo, times, parameters, b, d, J, T, Omega){
	
	D <- length(Xo) # dimension
	n_covs <- D+(D^2-D)/2 # Number of vars/covs: upper triangle + diagonal
	Mo <-  numeric(n_covs)	# initialize covariances at 0
	yo = c(Xo/Omega, Mo)
	eqns <- function(t,y,p){ linnoise(t,y,p, b, d, J, T) }

	out <- lsoda(yo, times, eqns, parameters)
	len <- dim(out)[2]
#	out[ , (D+2):(dim(out)[2]) ] <- out[ , (D+2):len ]*Omega
	out[,2:len] <- Omega*out[,2:len]
	out
}



# Transition=function(){0}
meanonly <- function(t,y,p, birth, death, Jacobian, Transition){
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
	net_transitions <- sapply(1:D, 
		function(i){
			sum(transition_matrix[,i]-transition_matrix[i,]) 
		})
	
	# macroscopics
	yd[1:D] <- birth(t,y,p) - death(t,y,p)  + net_transitions

	list(yd)
}

mean_only <- function(Xo, times, parameters, b, d, J, T, Omega){
	
	D <- length(Xo) # dimension
	n_covs <- D+(D^2-D)/2 # Number of vars/covs: upper triangle + diagonal
	Mo <-  numeric(n_covs)	# initialize covariances at 0
	yo = c(Xo/Omega, Mo)
	eqns <- function(t,y,p){ meanonly(t,y,p, b, d, J, T) }

	out <- lsoda(yo, times, eqns, parameters)
	len <- dim(out)[2]
#	out[ , (D+2):(dim(out)[2]) ] <- out[ , (D+2):len ]*Omega
	out[,2:len] <- Omega*out[,2:len]
	out
}







#depricated way
old_linear_noise_approx <- function(Xo, times, parameters, b, d, J, T, Omega){
	
	D <- length(Xo) # dimension
	Mo <-  numeric(D+(D^2-D)/2)	# initialize correct size covariance matrix
	yo = c(Xo/Omega, Mo)

	# extra function call slows this down, but keeps logic intuitive
	Nb <- function(t,y,p){b(t,y,p)/Omega }
	Nd <- function(t,y,p){d(t,y,p)/Omega }
	NJ <- function(t,y,p){J(t,y,p)/Omega }
	NT <- function(t,y,p){T(t,y,p)/Omega }
	eqns <- function(t,y,p){ linnoise(t,y,p, Nb, Nd, NJ, NT) }

	out <- lsoda(yo, times, eqns, parameters)
	print(mean(out[,4]))
	out[ , 2:(D+1) ] <- Omega*out[ , 2:(D+1) ]
	out[ , (D+2):(dim(out)[2]) ] <- Omega^2*out[ , (D+2):(dim(out)[2]) ]
	print(mean(out[,4]))
	out
}


