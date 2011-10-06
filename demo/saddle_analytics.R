# saddle_analytics.R
params <- c(K=1000, e=.5, h=200, p=2)
b <- function(x, params){ params["e"]*params["K"]*x^params["p"]/(x^params["p"]+params["h"]^params["p"]) }

## "Changing Mean" Dynamics
d <- function(x, a, params){ params["e"]*x+a}
A<- seq(50,250, length=100)


x <- seq(0, 1300, by=5)

require(rootSolve)
xhat <- function(a, params, interval=c(0,1000)){
		F <- function(x) {b(x, params) - d(x, a, params) }
		uniroot.all(F, interval=interval)
		
	}
	
bifur <- function(){
	for(a in A){
		plot(x, b(x,params)-d(x,a, params), type="l", ylim=c(-150, 900), ylab="rate")
		lines(x, b(x,params), lwd=4, lty=1, col="darkblue")
		lines(x, d(x,a, params), lwd=4, lty=1, col="darkred")
		abline(h=0, lwd=4, lty=2, col="darkgrey")
		crit_pts <- sort(xhat(a, params))
#		print(crit_pts)
		if(length(crit_pts) == 1){
			points(crit_pts, 0, pch=19, cex=1.5)
		} else if(length(crit_pts > 1)){
			points(min(crit_pts), 0, cex=1.5)
			points(max(crit_pts), 0, pch=19, cex=1.5)
		}
	}
}
require(animation)
## GIF animation, loop=0 for infinite loop
#saveMovie(bifur(), interval=.2, loop=0, outdir=getwd())

## shockwave flash animation:
saveSWF(bifur(), interval=.2, dev="png", swfname="bifur.swf", height=300, width=300)

lambda <- function(x, a, params){ params["p"]*params["e"]*params["K"]*x^(params["p"]-1)/(x^params["p"]+params["h"]^params["p"]) - params["p"]*params["e"]*params["K"]*x^(2*params["p"]-1)/(x^params["p"]+params["h"]^params["p"])^2-params["e"] }


crit_pts <- sapply(A, function(a){ xhat(a, params) })
good_xhat <-sapply(1:length(crit_pts), function(i) max(crit_pts[[i]]) )
alpha <- sapply(1:length(good_xhat), function(i) lambda(good_xhat[i], A[i], params))
sigma_squared <- sapply(1:length(good_xhat), function(i) b(good_xhat[i], params)+d(good_xhat[i], A[i], params))

png("model1.png", 1000, 1000)
par(mfrow=c(2,2))
plot(A, good_xhat, type="l", lwd=4,xlab="bifurcation parameter a", ylab="x hat", cex.lab=2, cex.axis=2 )
plot(A, alpha, type="l", lwd=4, xlab="bifurcation parameter a", cex.lab=2, cex.axis=2 )
plot(A, sigma_squared, type="l", lwd=4, xlab="bifurcation parameter a", cex.lab=2, cex.axis=2 )
plot(A, -sigma_squared/(2*alpha), type="l", lwd=4, xlab="bifurcation parameter a", cex.lab=2, cex.axis=2 )
dev.off()



## "Constant Mean" Dynamics
d <- function(x, a, params) { a*(x-params["K"])+params["e"]*params["K"] }
A<- seq(1,0,length=100)
saveSWF(bifur(), interval=.2, dev="png", swfname="bifur_mean_const.swf", height=300, width=300)


lambda <- function(x, a, params){ params["p"]*params["e"]*params["K"]*x^(params["p"]-1)/(x^params["p"]+params["h"]^params["p"]) - params["p"]*params["e"]*params["K"]*x^(2*params["p"]-1)/(x^params["p"]+params["h"]^params["p"])^2-a }

crit_pts <- sapply(A, function(a){ xhat(a, params) })
good_xhat <-sapply(1:length(crit_pts), function(i) max(crit_pts[[i]]) )
alpha <- sapply(1:length(good_xhat), function(i) lambda(good_xhat[i], A[i], params))
sigma_squared <- sapply(1:length(good_xhat), function(i) b(good_xhat[i], params)+d(good_xhat[i], A[i], params))



png("model2.png", 1000, 1000)
par(mfrow=c(2,2))
plot(A, good_xhat, type="l", lwd=4,xlab="bifurcation parameter a", ylab="x hat", cex.lab=2, cex.axis=2  )
plot(A, alpha, type="l", lwd=4, xlab="bifurcation parameter a", cex.lab=2, cex.axis=2 )
plot(A, sigma_squared, type="l", lwd=4, xlab="bifurcation parameter a", cex.lab=2, cex.axis=2 )
plot(A, -sigma_squared/(2*alpha), type="l", lwd=4, xlab="bifurcation parameter a", cex.lab=2, cex.axis=2 )
dev.off()




## Note that updating function definition updates prior functions, so the above doesn't need to redefine bifur():
f <- function(x){ x }
g <- function(x) { 2*f(x) }
f <- function(x) {2*x }
g(2)






## Uniroot.all makes the previous code here unnecssary 
