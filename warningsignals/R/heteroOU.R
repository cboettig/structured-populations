
heteroOU <- function(X, pars){
	if(is(X, "ts")){
		t <- time(X)
		x <- X@.Data
		N<- length(X)
	} else {
		t <- X[,1]
		x <- X[,2]
		N <- length(x)
	}
## Rescaling 
#	T <- max(t)
#	t <- t/T
#	pars[1] <- pars[1]*T
#	pars[2] <- pars[2]*T
#	pars[4] <- pars[4]*T
	lik <- 0.0
	out <- .C("heteroOU", as.double(lik), as.double(pars), as.double(x), as.double(t), as.integer(N))

	if(pars["Ro"] < 0){
		out[[1]] <- -Inf 
	}
	if(pars["sigma"] < 0){
		print("WHoooa!")
		out[[1]] <- -Inf
	}
	-out[[1]]
}

#pars <- c(0,0,0,1)
#X <- matrix(c(1:100, rnorm(100)), ncol=2)
#X <- ts(rnorm(1000))
#heteroOU(X,pars)
