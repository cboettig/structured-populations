# noam.R
n <- 20
sigma<- 2
alpha<-.05
X <- sigma*rnorm(n)
thetas <- rep(0, n)		# stable state of each patch
sigmas <- rep(sigma,n)
d<- .5
dt<- .05


update <-function(X, alpha){	
	X +
	(d*c(X[-1], X[1]) +d*c(X[n], X[-n]))*dt  +	# Diffusion (on S1)
	alpha*(thetas - X)*dt +		# return process
	sigmas*rnorm(n)*sqrt(dt) 		# noise
}


data <- sapply(1:500, function(i){
	alpha <- 5-.01*i
	X<<- update(X, alpha)
})

normalize_data <- function(data){
	sapply(1:dim(data)[2], function(i) data[,i]/sd(data[,i]) )
}

require(plotrix)
#color2D.matplot(data, border=NA)
color2D.matplot(normalize_data(data), border=NA)

