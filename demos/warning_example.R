# Example using data simulated from the models themselves
reps <- 50

require(stochPop)
ou = list(alpha=2, theta=3, sigma=1)
wa = list(alpha_0=2, theta=3, sigma=1, beta=2)
cpt = list(alpha_1=2, alpha_2=.2, theta=3, sigma=1, t_shift=5)
class(ou) <- c("OU", "list")
class(wa) <- c("warning", "list")
class(cpt) <- c("changePt", "list")
Dt <- 1
Xo <- 3


# Simulate an example dataset without warning signal 
ou$data <- simulate(pars=ou, T = 10, N=2000, X0=Xo)
wa$data <- simulate.warning(pars=wa, T = 10, N=2000, X0=Xo)
cpt$data <- simulate.changePt(pars=cpt, T = 10, N=2000, X0=Xo)


# Fit the models, even MLE call works happily with Nelder-Mead algorithm
ou <- update.OU(pars = ou, X = ou$data, method="L") 
wa <- update.warning(wa, X = wa$data)
cpt <- update.changePt(cpt, X = cpt$data)


ou_boot <- bootstrap(ou, reps = reps)
wa_boot <- bootstrap(wa, reps = reps)
cpt_boot <- bootstrap(cpt, reps = reps)

save(list=ls(), file= "warning_example.Rdat") 
png("ou_boot",800,200); plot_bootstrap(ou_boot); dev.off()
png("wa_boot",800,200); plot_bootstrap(wa_boot); dev.off()
png("cpt_boot",800,200); plot_bootstrap(cpt_boot); dev.off()










