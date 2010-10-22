# saddle_analytics.R
params <- c(K=1000, e=.5, h=200, p=2)
b <- function(x, params){ params["e"]*params["K"]*x^params["p"]/(x^params["p"]+params["h"]^params["p"]) }
d <- function(x, a, params){ params["e"]*x+a}

lambda <- function(x, a, params){ params["p"]*params["e"]*params["K"]*x^(params["p"]-1)/(x^params["p"]+params["h"]^params["p"]) - pparams["p"]*params["e"]*params["K"]*x^(2*params["p"]-1)/(x^params["p"]+params["h"]^params["p"])^2-params["e"] }


x <- seq(0, 1000, by=5)
a<- 100
#par(mfrow=c(2,1))
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





# The canonical form of the bifurcation is second order,
# dx/dt = x^2 + a
# and involves a symmetric approach.  
xhat_left <- function(a){ 
	f <- function(x){ x^2+a }
	uniroot(f, c(-100,0))$root
}
xhat_right <- function(a){ 
	f <- function(x){ x^2+a }
	uniroot(f, c(0, 100))$root
}

a <- seq(-5,0,length=100)
yl <- sapply(a, xhat_left)
yr <- sapply(a, xhat_right)
plot(c(a,a), c(yl,yr), pch=19, xlab="lambda", ylab="x_hat")


