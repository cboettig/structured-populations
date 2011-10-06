
## Example using the saddle-node simulation
reps <- 100
cpu <-2 


# Simulate a dataset from the full individual, nonlinear model
T<- 40*25
n_pts <- 201
require(stochPop)
pars = c(Xo = 730, e = 0.5, a = 100, K = 1000, h = 200, 
    i = 0, Da = .12, Dt = 0, p = 2)
sn <- saddle_node_ibm(pars, times=seq(0,T, length=n_pts))
X <- ts(sn$x1,start=sn$time[1], deltat=sn$time[2]-sn$time[1])

plot(X, cex.lab=2, cex.axis=2, lwd=3, xlab="time", ylab="pop")



Y <- loess.smooth(time(X), X, evaluation=length(X))
Y <- ts(Y$y, start=sn$time[1], deltat=sn$time[2]-sn$time[1])
Y <- X-Y
#plot(X-Y)
#X <- X-Y

# Initialize some parameter estimates as initial guesses for fitting 
theta <- mean(X)
sigma <- sd(X)
alpha <- sigma

# Intialize the three models we'll test
ou <- list(alpha=alpha, theta=theta, sigma=sigma)
wa <- list(alpha_0=alpha, theta=theta, sigma=sigma, beta=0)
cpt <- list(alpha_1=alpha, alpha_2=alpha, theta=theta, sigma=sigma, t_shift=time(X)[length(X)]/2)
class(ou) <- c("OU", "list")
class(wa) <- c("warning", "list")
class(cpt) <- c("changePt", "list")

# Fit the models, consider MLE call for OU & warning
ou <- update(ou, X = X) 
ou <- update(ou, X = X) 
wa <- update(wa, X = X)
wa <- update(wa, X = X)
cpt <- update(cpt, X = X)
cpt <- update(cpt, X = X)

print(c(ou$loglik, wa$loglik, cpt$loglik))


#model_list <- list(ou=ou, wa=wa, cpt=cpt)
#model_boots <- bootstrapLR(model_list, rep=reps, cpu=cpu)


## Plot the results of the model choice
#LRplot(model_boots, test=2,null=1, main= "1 vs 2")


pt <- power_pair(wa, ou, nboot=reps, cpu=cpu)

save(list=ls(), file="warning4.Rdat")
## And the bootstraps of the parameter estimates
#png("ou_boot.png",1600,400); 
#plot_bootstrap(model_boots, model=1, cex.lab=2, cex.axis=2, lwd=3); 
#dev.off()
#png("wa_boot.png",1600,400); plot_bootstrap(model_boots, model=2, cex.lab=2, cex.axis=2, lwd=3); dev.off()
#png("cpt_boot.png",1600,400); plot_bootstrap(model_boots, model=3, cex.lab=2, cex.axis=2, lwd=3); dev.off()




#pt <- power_test(ou, wa, nboot=reps, cpu=cpu, values=c(.001, .01, .5, seq(1,15, by=3)) )
#save(list=ls(), file= "warning_example3.Rdat") 


#png("power.png")
#plot(pt$values, pt$power)
#dev.off()
