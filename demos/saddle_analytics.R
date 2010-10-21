# saddle_analytics.R
params <- c(K=1000, e=.5, h=200)
p <- 5
b <- function(x, params){ params["e"]*params["K"]*x^p/(x^p+params["h"]^p) }
d <- function(x, a, params){ params["e"]*x+a}

lambda <- function(x, a, params){ p*params["e"]*params["K"]*x^(p-1)/(x^p+params["h"]^p) - p*params["e"]*params["K"]*x^(2*p-1)/(x^p+params["h"]^p)^2-params["e"] }


x <- seq(0, 1000, by=5)
a<- 200
plot(x, b(x,params)-d(x,a, params), type="l", ylim=c(-150, 900))
lines(x, b(x,params), lwd=4, lty=1, col="darkblue")
lines(x, d(x,a, params), lwd=4, lty=1, col="darkred")


xhat <- function(a){
	F <- function(x, a, params){abs(b(x, params) - d(x, a, params)) } 
	nlm(F, 1000, a, params)$estimate
}

a <- seq(100, 300, length=500)
X <- sapply(a, xhat)

par(mar=c(5,4,4,5)+.1)
plot(a, X, xlab="a", ylab="Pop Equilibrium", pch=17)
par(new=TRUE)
plot(a, lambda(X, a, params), xlab="", ylab="", pch=19, xaxt="n", yaxt="n", col="red")
axis(4, col="red", col.ticks="red")
mtext("lambda", side=4, line=3, col="red")
legend("left",col=c("black","red"),pch=c(17,19),legend=c("Population","lambda"))
abline(h=0, col="red")






x <- 300:800
plot(x, b(x,params)-d(x,174,params) )


xhat_uni <- function(a){ 
	f <- function(x){ b(x, params) - d(x, a, params) }
	uniroot(f, c(200, 1000))$root
}
xhat_uni(100)
plot(a, xhat_uni(a))


