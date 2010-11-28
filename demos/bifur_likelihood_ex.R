# bifur_likelihood_ex.R

source("../R/likelihood_bifur_models.R")
pars = c(r=10, theta=3, beta=2)
m <- init_sdemodel(pars =pars, Xo = 6.2, model="SN", N=200)
X <- simulate.SN(m)
# set the initial search values (to be other than the true ones)
m$pars <- c(r=11,theta=1.5,beta=1)

out <- update.SN(m, X)
print(out$pars)


## Figure 
png("saddle_node_fit.png", width=400, height=400)
curve( -(x-theta)^2+r, 0, 10, ylim=c(-2, r+1), lwd=3, main="true vs estimated model", col="darkgray")
curve( -(x-out$pars["theta"])^2+out$pars["r"], 0, 10, ylim=c(-2, r+1),add=T, col="red", lty=2, lwd=3)
text(1, 2, paste("est: ", "r = ", as.character(round(out$pars[1],2)), "theta = ", as.character(round(out$pars[2],2)), "beta = ", as.character(round(out$pars[3],2))), pos=4)
text(1, 1, paste("true: ", "r = ", as.character(pars[1]), "theta = ", as.character(pars[2]), "beta = ", as.character(pars[3])), pos=4)
abline(h=0, lty=2)
par(new=TRUE)
plot(X@.Data, time(X),,type="l",col="blue",xlim=c(0,10), xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("data",side=4,line=3)
dev.off()

gitcom <- system('git commit -a -m "autocommit"', intern=TRUE)
# only first line of git commit commandline return will be printed in flickr discription
system(paste('flickr_upload --tag="stochpop bifurcation"', '--description="', gitcom, '"',  ' saddle_node_fit.png'))
system(paste('hpc-autotweets ', gitcom, ' done'))


A <- seq(1,10, by=.1)
l <- sapply(A, function(a){	
	pars['theta'] <- a
	SN.lik(X, pars)
})
png("likelihood_cross-section.png")
plot(A, l, xlab="theta", ylab="-loglik")
lines(A,l)
dev.off()
system('flickr_upload --tag="stochpop likelihood saddle-node" likelihood_cross-section.png')

