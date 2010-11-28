# bifur_likelihood_ex.R

source("../R/likelihood_bifur_models.R")
r<-10
theta<-3
beta<-1
m <- init_sdemodel(pars = c(r=r, theta=theta, beta=beta), Xo = 6.2, model="SN", N=200)
X <- simulate.sdemodel(m)
# IC
m$pars <- c(r=15,theta=5,beta=.1)

out <- update.sn(m, X)
print(out$pars)

png("saddle_node_fit.png", width=800, height=400)
par(mfrow=c(1,2))
curve( -(x-theta)^2+r, 0, 10, ylim=c(-2, r+1), lwd=3, main="true vs estimated model")
curve( -(x-out$pars["theta"])^2+out$pars["r"], 0, 10, ylim=c(-2, r+1),add=T, col="red", lty=2, lwd=3)
plot(X, lwd=3, main="data")
dev.off()
