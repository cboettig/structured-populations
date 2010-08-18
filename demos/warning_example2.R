
## Example using the saddle-node simulation
reps <- 100
cpu <- 2

T<- 200
require(stochPop)
pars = c(Xo = 570, e = 0.5, a = 100, K = 1000, h = 200, 
    i = 0, Da = .5, Dt = 0)
sn <- saddle_node_ibm(pars, times=seq(0,T, length=2000))
X <- ts(sn$x1,start=0, end=T )
plot(X)

theta <- mean(X)
sigma <- sd(X)
alpha <- sigma

ou <- list(alpha=alpha, theta=theta, sigma=sigma)
wa <- list(alpha_0=alpha, theta=theta, sigma=sigma, beta=2)
cpt <- list(alpha_1=alpha, alpha_2=alpha, theta=theta, sigma=sigma, t_shift=time(X)[length(X)]/2)
class(ou) <- c("OU", "list")
class(wa) <- c("warning", "list")
class(cpt) <- c("changePt", "list")
Dt <- 1
Xo <- 3

# Fit the models, consider MLE call for OU & warning
ou <- update.OU(pars = ou, X = X) 
wa <- update.warning(wa, X = X)
cpt <- update.changePt(cpt, X = X)
model_list <- list(ou=ou, wa=wa, cpt=cpt)

model_boots <- bootstrapLR(model_list, rep=reps, cpu=cpu)


#cpt_boot <- bootstrap(cpt, reps = reps, cpu=cpu)
save(list=ls(), file= "warning_example.Rdat") 

png("model_choice.png", 1600, 400)
par(mfrow=c(1,3))
LRplot(model_boots, 1,2, main= "OU vs warning")
LRplot(model_boots, 1,3, main="OU vs ChangePt")
LRplot(model_boots, 2,3, main = "warning vs ChangePt" )
dev.off()


png("ou_boot",1600,400); plot_bootstrap(model_boots, model=1); dev.off()
png("wa_boot",1600,400); plot_bootstrap(model_boots, model=2); dev.off()
png("cpt_boot",1600,400); plot_bootstrap(model_boots, model=3); dev.off()


