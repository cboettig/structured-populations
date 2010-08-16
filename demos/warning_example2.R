
## Example using the saddle-node simulation
reps <- 50
cpu <- 2

T<- 40
require(stochPop)
pars = c(Xo = 570, e = 0.5, a = 160, K = 1000, h = 200, 
    i = 0, Da = 0.0, Dt = 100)
sn <- saddle_node_ibm(pars, times=seq(0,T, length=40))
X <- ts(sn$x1,start=0, end=T )
plot(X)

theta <- mean(X)
sigma <- sd(X)
alpha <- sigma

ou = list(alpha=alpha, theta=theta, sigma=sigma)
wa = list(alpha_0=alpha, theta=theta, sigma=sigma, beta=2)
cpt = list(alpha_1=alpha, alpha_2=alpha, theta=theta, sigma=sigma, t_shift=time(X)/2)
class(ou) <- c("OU", "list")
class(wa) <- c("warning", "list")
class(cpt) <- c("changePt", "list")
Dt <- 1
Xo <- 3

# Fit the models, even MLE call works happily with Nelder-Mead algorithm
ou <- update.OU(pars = ou, X = X) 
wa <- update.warning(wa, X = X)
#cpt <- update.changePt(cpt, X = X)


ou_boot <- bootstrap(ou, reps = reps, cpu=cpu)
wa_boot <- bootstrap(wa, reps = reps, cpu=cpu)
cpt_boot <- bootstrap(cpt, reps = reps, cpu=cpu)

save(list=ls(), file= "warning_example.Rdat") 

png("ou_boot",800,200); plot_bootstrap(ou_boot); dev.off()
png("wa_boot",800,200); plot_bootstrap(wa_boot); dev.off()
png("cpt_boot",800,200); plot_bootstrap(cpt_boot); dev.off()


