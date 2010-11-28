# bifur_likelihood_ex.R

source("../R/likelihood_bifur_models.R")
pars = c(r=5, theta=3, beta=1)
m <- init_sdemodel(pars =pars, Xo = 8, model="SN", N=500, T=100)
X <- simulate.SN(m)
# set the initial search values (to be other than the true ones)
m$pars <- c(r=21,theta=.5,beta=.1)

out <- update.SN(m, X)
print(out$pars)




## Figure 
png("saddle_node_fit.png")
curve( -(x-pars['theta'])^2+pars['r'], 0, 10, ylim=c(-2, r+1), lwd=3, main="true vs estimated model", col="darkgray")
curve( -(x-out$pars["theta"])^2+out$pars["r"], 0, 10, ylim=c(-2, r+1),add=T, col="red", lty=2, lwd=3)
#text(1,4, paste("N = ", m$N, " T = ", m$T),pos=4)
text(1, 3, paste("est: ", "r = ", as.character(round(out$pars[1],2)), "theta = ", as.character(round(out$pars[2],2)), "beta = ", as.character(round(out$pars[3],2))), pos=4)
text(1, 2, paste("init: ", "r = ", as.character(m$pars[1]), "theta = ", as.character(m$pars[2]), "beta = ", as.character(m$pars[3])), pos=4)
text(1, 1, paste("true: ", "r = ", as.character(pars[1]), "theta = ", as.character(pars[2]), "beta = ", as.character(pars[3])), pos=4)
abline(h=0, lty=2)
par(new=TRUE)
plot(X@.Data, time(X),,type="l",col="blue",xlim=c(0,10), xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("data",side=4,line=3)
dev.off()

# only first line of git commit will be used 
gitcom <- system('git commit -a -m "autocommit"', intern=TRUE)
system(paste('flickr_upload --tag="stochpop bifurcation"', '--description="', gitcom, '"',  ' saddle_node_fit.png'))
system(paste('hpc-autotweets "', gitcom, ' done"'))





ou_pars <- list(alpha=1, theta=1, sigma=1)
ou <- update.OU(ou_pars, X)
r <- ou$alpha^2/4
theta <- ou$theta - ou$alpha/2
beta <- ou$sigma^2/(2*ou$alpha*sqrt(r))



linearmodel <- m
linearmodel$pars <- c(r=r, theta=theta, beta=beta)
print(SN.lik(X, linearmodel))
linearfit <- update.SN(linearmodel, X)

## Plot the likelihood surface 
A <- seq(1,10, length=200)
l <- sapply(A, function(a){	
	pars['theta'] <- a
	SN.lik(X, pars)})
png("likelihood_cross-section.png")
plot(A, l, xlab="theta", ylab="-loglik")
lines(A,l)
dev.off()
system('flickr_upload --tag="stochpop likelihood saddle-node" likelihood_cross-section.png')



