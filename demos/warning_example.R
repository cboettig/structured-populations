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


# Fit the models, even MLE call works happily with Nelder-Mead algorithm
fit <- OU.fitML(X1, pars_ou) 
fit2 <- warning.fitML(X2, pars_wa)
fit3 <- changePt.fitML(X3, pars_cpt)



bootstrap <- function(fit_model, reps=100, cpu=2){
	require(snowfall)
	if(cpu>1){
		sfInit()
	} else {
		sfinit(parallel=TRUE, cpu=cpu)
		sfLibrary(stochPop)
		sfExportAll()
	}

	data <- sfSapply(1:reps, function(i) )

}

## Example using the saddle-node simulation
sn <- saddle_node_ibm()
X <- ts(sn$x1)













