require(stochPop)
ou = list(alpha=2, theta=3, sigma=1)
wa = list(alpha_0=2, theta=3, sigma=1, beta=2)
cpt = list(alpha_1=2, alpha_2=.2, theta=3, sigma=1, t_shift=5)
class(ou) <- c("OU", "list")
class(wa) <- c("warning", "list")
class(cpt) <- c("changePt", "list")
Dt <- 1
Xo <- 3

source("R/sde_likelihood.R")

# Simulate an example dataset without warning signal 
pars_ou$data <- simulate(pars=ou, T = 10, N=2000, X0=Xo)
pars_wa$data <- simulate.warning(pars=wa, T = 10, N=2000, X0=Xo)
pars_cpt$data <- simulate.changePt(pars=cpt, T = 10, N=2000, X0=Xo)


# Fit the models, even MLE call works happily with Nelder-Mead algorithm
ou <- update.OU(pars = ou, X = ou$data, method="L") 
wa <- update.warning(wa, X = wa$data)
cpt <- update.changePt(cpt, X = cpt$data)


## Example using the saddle-node simulation
sn <- saddle_node_ibm()
X <- ts(sn$x1)













