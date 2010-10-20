# saddle_analytics.R
params <- c(K=1000, e=.5, h=200)

b <- function(x, params){ params["e"]*params["K"]*x^2/(x^2+params["h"]^2) }
d <- function(x, a, params){ params["e"]*x+a}

lambda <- function(x, a, params){ 2*params["e"]*params["K"]*x/(x^2+h^2) - 2*params["e"]*params["K"]*x^3/(x^2+h^2)^2-e }


x <- seq(0, 1000, by=5)
plot(x, b(x,params)-d(x,100, params), type="l")

a <- seq(100, 300, by=10)

xhat <- function(a){
	F <- function(x, a, params){abs(b(x, params) - d(x, a, params)) } 
	nlm(F, 1000, a, params)$estimate
}

plot(a, sapply(a, xhat))



f <- function(x){ b(x, params) - d(x, a, params) }
fdash <- function(x){ lambda(x, 100, params) }

require(elliptic)
newton.rapheson(1000, f, dash)

xhat <- function(a){ 
	f <- function(x){ b(x, params) - d(x, a, params) }
	uniroot(f, c(200, 1000))$root
}

plot(a, xhat(a))


