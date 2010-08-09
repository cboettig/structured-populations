require(stochPop)
set.seed(123)
theta <- 3
alpha <- 1
sigma <- 2
alpha_0 <- 1
beta <- .2
pars = list(alpha_0=alpha, theta=theta, sigma=sigma, beta=beta)

Dt <- 1
Xo <- 1

# Numeric and analytic approximations to the variance agree
numeric_V(Dt, pars)
analytic_V(Dt, pars)

# Compare the probability density between the two models.  Should roughly agree with small beta,
# and should agree exactly with very small beta.  
warning_model(Dt, Xo, pars)
setOU(Dt, Xo, theta=c(theta*alpha, alpha, sigma) )

# Simulate an example dataset without warning signal 
X <- sde.sim(model="OU", theta= c(theta*alpha,alpha,sigma), X0=Xo, N=1000, T=10) # can specify N & delta and will calc T, or T & N.  

t_shift <- max(time(X))/2

# Evaluate the likelihood of the models with true parameter values
OU.lik(alpha*theta, alpha, sigma) 
warning.lik(alpha_0, theta, sigma, beta)
TwoRates.lik(alpha, alpha, theta, sigma, 5)


# Fit the models

# Parametrization dXt = (O1 - O2*Xt)dt + O3*dWt
# alpha = theta2, theta = theta1/theta2, sigma = theta3

fit <- mle(OU.lik, start=list(theta1=1, theta2=.5, theta3=.5), method="L-BFGS-B", lower=c(-Inf, 0,0))
summary(fit)
fit2 <- mle(warning.lik, start=list(alpha_0=.5, theta=1, sigma=5, beta=2), method="L-BFGS-B", lower=c(0,0,0,1e-9))
summary(fit2)
fit3 <- mle(TwoRates.lik, start=list(alpha_1=.5, alpha_2 = .5, theta = 1, sigma =3, t_shift = 5), method="L-BFGS-B", lower=c(0,0,-Inf, 0, 0), upper=c(Inf,Inf, Inf,Inf, max(time(X))-2*deltat(X) ))





## formatted directly for optim call, note data is always assumed to be in X
OU.likfn <- function(x){
		OU.lik(x[1], x[2], x[3])
}
TwoRates.likfn <- function(x){
		TwoRates.lik(x[1], x[2], x[3], x[4], x[5])
}


o <- optim(c(alpha*theta, alpha, sigma), OU.likfn)
o3 <- optim(c(alpha, alpha, theta, sigma, t_shift), TwoRates.likfn, method="L", lower=c(0,0,0,0,0))


