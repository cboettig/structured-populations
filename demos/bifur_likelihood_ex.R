# bifur_likelihood_ex.R
## Figure plotting function
source("../R/likelihood_bifur_models.R")

plotcurves <- function(m, out, pars, X, oucurve=NULL){ # m is the initial conditions model (input to update.SN), out the fit model (output of update.SN), pars the parameters used for the simulated data, X
	curve( -(x-pars['theta'])^2+pars['r'], 0, 10, ylim=c(-2, r+1), lwd=3, main="true vs estimated model", col="darkgray")
	curve( -(x-out$pars["theta"])^2+out$pars["r"], 0, 10,add=T, col="red", lty=2, lwd=3)
	if(!is.null(oucurve)){
		curve( -(x-oucurve$pars["theta"])^2+oucurve$pars["r"], 0, 10,add=T, col="green", lty=2, lwd=3)
	}
#text(1,4, paste("N = ", m$N, " T = ", m$T),pos=4)
	text(1, 3, paste("est: ", "r = ", as.character(round(out$pars[1],2)), "theta = ", as.character(round(out$pars[2],2)), "beta = ", as.character(round(out$pars[3],2))), pos=4)
	text(1, 2, paste("init: ", "r = ", as.character(m$pars[1]), "theta = ", as.character(m$pars[2]), "beta = ", as.character(m$pars[3])), pos=4)
	text(1, 1, paste("true: ", "r = ", as.character(pars[1]), "theta = ", as.character(pars[2]), "beta = ", as.character(pars[3])), pos=4)
	abline(h=0, lty=2)
	par(new=TRUE)
	plot(X@.Data, time(X),,type="l",col=rgb(0,0,1,.4),xlim=c(0,10), xaxt="n",yaxt="n",xlab="",ylab="")
	axis(4)
	mtext("data",side=4,line=3)
}


pars = c(r=5, theta=3, beta=.1)
m <- init_sdemodel(pars =pars, Xo = 8, model="SN", N=500, T=100)
X <- simulate.SN(m)
# set the initial search values (to be other than the true ones)
m$pars <- c(r=21,theta=.5,beta=.01)

out <- update.SN(m, X)


#plot 
png("saddle_node_fit.png")
plotcurves(m, out, pars, X)
dev.off()

# only first line of git commit will be used 
gitcom <- system('git commit -a -m "autocommit"', intern=TRUE)
system(paste('flickr_upload --tag="stochpop bifurcation"', '--description="', gitcom, '"',  ' saddle_node_fit.png'))
system(paste('hpc-autotweets "', gitcom, ' done"'))



## estimate the linear model directly and convert parameters to the quadratic SN model
ou_pars <- list(alpha=1, theta=1, sigma=1)
ou <- update.OU(ou_pars, X)
r <- ou$alpha^2/4
theta <- ou$theta - ou$alpha/2
beta <- ou$sigma^2/(2*ou$alpha*sqrt(r))
linearmodel <- m
linearmodel$pars <- c(r=r, theta=theta, beta=beta)

png("linear_fit.png")
plotcurves(m=m, out=out, oucurve=linearmodel, pars=pars, X=X)
legend("topright", c("true", "quad est", "transform lin est"), col=c("darkgray", "red", "green"), lty=1 )
dev.off()
gitcom <- system('git commit -a -m "autocommit"', intern=TRUE)
system(paste('flickr_upload --tag="stochpop bifurcation"', '--description="', gitcom, '"',  ' linear_fit.png'))


## fit the quadratic model using the linear estimate as a starting point
#linearfit <- update.SN(linearmodel, X)





## Plot the likelihood surface 
#A <- seq(1,10, length=200)
#l <- sapply(A, function(a){	
#	pars['theta'] <- a
#	SN.lik(X, pars)})
#png("likelihood_cross-section.png")
#plot(A, l, xlab="theta", ylab="-loglik")
#lines(A,l)
#dev.off()
#system('flickr_upload --tag="stochpop likelihood saddle-node" likelihood_cross-section.png')




