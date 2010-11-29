# bifur_likelihood_ex.R
source("../R/likelihood_bifur_models.R")
source("../R/sde_likelihood.R")

## Set up a quadratic model and simulate a dataset 
pars = c(r=5, theta=3, beta=.8)
m <- init_sdemodel(pars =pars, Xo = 8, model="SN", N=500, T=100)
X <- simulate.SN(m)

## estimate the linear model directly and convert parameters to the quadratic SN model
ou_pars <- list(alpha=1, theta=1, sigma=1)
ou <- update.OU(ou_pars, X)
r <- ou$alpha^2/4
theta <- ou$theta - ou$alpha/2
beta <- ou$sigma^2/(2*ou$alpha*sqrt(r))

## use the linear model as initial conditions for quadratic fit  
m$pars <- c(r=r, theta=theta, beta=beta)
linearmodel <- m


## perform the quadratic fit
out <- update.SN(m, X)

#plot 
png("saddle_node_fit.png")
plotcurves(m=m, out=out, oucurve=linearmodel, pars=pars, X=X)
legend("topright", c("true", "quad est", "transform lin est"), col=c("darkgray", "red", "green"), lty=1 )
dev.off()

##social reporting
gitcom <- system('git commit -a -m "autocommit"', intern=TRUE)[[1]]
system(paste('flickr_upload --tag="stochpop bifurcation"', '--description="', gitcom, '"',  ' saddle_node_fit.png'))
system(paste('hpc-autotweets " run ', gitcom, 'done "'))


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




