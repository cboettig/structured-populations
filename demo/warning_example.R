# Example using data simulated from the models themselves
reps <- 200
cpu <- 2

require(stochPop)
ou = list(alpha=30, theta=300, sigma=30, T=10, N=100, X0=300)
wa = list(alpha_0=30, theta=300, sigma=30, beta=0, T=10, N=100, X0=300)
cpt = list(alpha_1=50, alpha_2=5, theta=300, sigma=30, t_shift=5, T=10, N=100, X0=300)
class(ou) <- c("OU", "list")
class(wa) <- c("warning", "list")
class(cpt) <- c("changePt", "list")
Dt <- 1

X <- simulate(cpt)
plot(X)
OU.lik(X, ou)
warning.lik(X, wa)



# Fit the models, consider MLE call for OU & warning
ou <- update.OU(pars = ou, X = X) 
wa <- update.warning(wa, X = X)
cpt <- update.changePt(cpt, X = X)
model_list <- list(ou=ou, wa=wa, cpt=cpt)

model_boots <- bootstrapLR(model_list, rep=reps, cpu=cpu)

save(list=ls(), file= "warning_example.Rdat") 

png("model_choice.png", 1600, 400)
par(mfrow=c(2,3))
LRplot(model_boots, test=2,null=1, main= "1 vs 2")
LRplot(model_boots, test=3, null=1, main="1 vs 3")
LRplot(model_boots, 3,2, main = "2 vs 3" )
LRplot(model_boots, 1,2)
LRplot(model_boots, 1,3)
LRplot(model_boots, 2,3)
dev.off()


png("ou_boot.png",1600,400); plot_bootstrap(model_boots, model=1, cex.lab=2, cex.axis=2, lwd=3); dev.off()
png("wa_boot.png",1600,400); plot_bootstrap(model_boots, model=2, cex.lab=2, cex.axis=2, lwd=3); dev.off()
png("cpt_boot.png",1600,400); plot_bootstrap(model_boots, model=3, cex.lab=2, cex.axis=2, lwd=3); dev.off()




