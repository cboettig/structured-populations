require(beetles)
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
X <- sde.sim(model="OU", theta= c(theta*alpha,alpha,sigma), X0=Xo, N=1000, T=1000) # can specify N & delta and will calc T, or T & N.  

# Evaluate the likelihood of the models with true parameter values
OU.lik(alpha*theta, alpha, sigma) 
warning.lik(alpha_0, theta, sigma, beta)


# Fit the models
fit <- mle(OU.lik, start=list(theta1=1, theta2=.5, theta3=.5), method="L-BFGS-B", lower=c(-Inf, 0,0))
summary(fit)
fit2 <- mle(warning.lik, start=list(alpha_0=.5, theta=1, sigma=2, beta=2), method="L-BFGS-B", lower=c(0,0,0,1e-9))
summary(fit2)


