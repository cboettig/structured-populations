# saddle_analytics.R
params <- c(K=1000, e=.5, h=200, p=2)
b <- function(x, params){ params["e"]*params["K"]*x^params["p"]/(x^params["p"]+params["h"]^params["p"]) }

## "Changing Mean" Dynamics
d <- function(x, a, params){ params["e"]*x+a}
A<- seq(50,250, length=100)

## "Constant Mean" Dynamics
#d <- function(x, a, params) { a*(x-params["K"])+params["e"]*params["K"] }
#A<- seq(1,0,length=100)

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
		print(crit_pts)
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
saveSWF(bifur(), interval=.2, dev="png", swfname="bifur.swf", height=300, width=300, outdir=getwd())




# "Changing mean dynamics"

lambda <- function(x, a, params){ params["p"]*params["e"]*params["K"]*x^(params["p"]-1)/(x^params["p"]+params["h"]^params["p"]) - pparams["p"]*params["e"]*params["K"]*x^(2*params["p"]-1)/(x^params["p"]+params["h"]^params["p"])^2-params["e"] }


## Uniroot.all makes the previous code here unnecssary 
