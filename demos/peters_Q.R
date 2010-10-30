#peters_Q.R

sim <- function(){
	N<-1000
	dX <- rnorm(N, 0, .01)
	X <- numeric(N)
	X[1] <- 0
	dY<- numeric(N)
	Y<- numeric(N)
	Y[1] <- 0
	for(i in 2:N) {
		X[i] <- X[i-1]+dX[i] 
		dY[i] <- rnorm(1, 0, abs(X[i]))
		Y[i] <- Y[i-1]+dY[i] 
	}
	Y[N]
}



data <- sapply(1:2000, function(i) sim())

(data)

plot(pdf(
