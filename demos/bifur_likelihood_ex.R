# bifur_likelihood_ex.R

source("../R/likelihood_bifur_models.R")
r<-10
theta<-3
beta<-1
curve( -(x-theta)^2+r, 0, 10, ylim=c(-2, r+1))
m <- init_sdemodel(pars = c(r=r, theta=theta, beta=beta), Xo = 10, model="SN")
X <- simulate.sdemodel(m)
plot(X)

out <- update.sn(m, X)
