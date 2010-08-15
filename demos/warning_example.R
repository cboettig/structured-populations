require(stochPop)
pars_ou = list(alpha=2, theta=3, sigma=1)
pars_wa = list(alpha_0=2, theta=3, sigma=1, beta=2)
pars_cpt = list(alpha_1=2, alpha_2=.2, theta=3, sigma=1, t_shift=5)
Dt <- 1
Xo <- 3

source("R/sde_likelihood.R")

# Simulate an example dataset without warning signal 
X1 <- OU.sim(T = 10, N=2000, X0=Xo, pars=pars_ou)
X2 <- warning.sim(T = 10, N=2000, X0=Xo, pars=pars_wa)
X3 <-changePt.sim(T = 10, N=2000, X0=Xo, pars=pars_cpt)


# Fit the models


fit <- OU.fitML(X1, pars_ou) 
fit2 <- warning.fitML(X2, pars_wa)
fit3 <- changePt.fitML(X3, pars_cpt)



## formatted directly for optim or mle call, which takes parameters to be optimized in named list as argument 




o <- optim(c(alpha*theta, alpha, sigma), OU.likfn)
bootstrap <- function(fit_model, reps=100, cpu=2){

}

## Example using the saddle-node simulation
sn <- saddle_node_ibm()
X <- ts(sn$x1)

oupars <- c( mean(X)*100, 100, sd(X) )
o <- optim( oupars, OU.likfn)
o2 <- optim( c(o$par[2], o$par[1]/o$par[2], o$par[3], 0 ), warning.likfn)












