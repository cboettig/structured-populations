require(stochPop)
pars = c(alpha_0=2, theta=3, sigma=1, beta=2)
pars2 = c(alpha_1=2, alpha_2=.2, theta=3, sigma=1, t_shift=0.5)
Dt <- 1
Xo <- 3

# Simulate an example dataset without warning signal 
X1 <- sde.sim(model="OU", theta= c(theta*alpha,alpha,sigma), X0=Xo, N=2000, T=10) # can specify N & delta and will calc T, or T & N.  
X2 <- warning.sim(T = 10, N=2000, X0=Xo, pars=pars)

X<-ts(X3)

# Evaluate the likelihood of the models with true parameter values
OU.lik(alpha*theta, alpha, sigma) 
warning.lik(alpha_0, theta, sigma, beta)
TwoRates.lik(alpha, alpha_2, theta, sigma, t_shift)


# Fit the models

# Parametrization dXt = (O1 - O2*Xt)dt + O3*dWt
# alpha = theta2, theta = theta1/theta2, sigma = theta3

fit <- mle(OU.lik, start=list(theta1=1, theta2=.5, theta3=.5), method="L-BFGS-B", lower=c(-Inf, 0,0))
summary(fit)
fit2 <- mle(warning.lik, start=list(alpha_0=.5, theta=1, sigma=5, beta=2), method="L-BFGS-B", lower=c(0,0,0,1e-9))
summary(fit2)
#fit3 <- mle(TwoRates.lik, start=list(alpha_1=.5, alpha_2 = .5, theta = 1, sigma =3, t_shift = 5), method="L-BFGS-B", lower=c(0,0,-Inf, 0, 0), upper=c(Inf,Inf, Inf,Inf, max(time(X))-2*deltat(X), fixed=list(alpha_2, start_t) ))





## formatted directly for optim call, note data is always assumed to be in X, little x just gives parameters
## note these functions all actually return the negative log likelihood
OU.likfn <- function(x){
		OU.lik(x[1], x[2], x[3])
}
warning.likfn <- function(x){
		warning.lik(x[1], x[2], x[3], x[4])
}
TwoRates.likfn <- function(x){
		TwoRates.lik(x[1], x[2], x[3], x[4], x[5])
}


o <- optim(c(alpha*theta, alpha, sigma), OU.likfn)
o3 <- optim(c(alpha, alpha, theta, sigma, t_shift), TwoRates.likfn, method="L", lower=c(0,0,0,0,0))

fit_ou
fit _warning
fit_TwoRates

bootstrap <- function(fit_model, reps=100, cpu=2){

}

## Example using the saddle-node simulation
sn <- saddle_node_ibm()
X <- ts(sn$x1)

oupars <- c( mean(X)*100, 100, sd(X) )
o <- optim( oupars, OU.likfn)
o2 <- optim( c(o$par[2], o$par[1]/o$par[2], o$par[3], 0 ), warning.likfn)
o3 <- optim( c(o$par[2], o$par[2], o$par[1]/o$par[2], o$par[3], max(time(X))/2 ), TwoRates.likfn)












