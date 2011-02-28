cpu <- 16
nboot <- 160

# CaCO3 record
caco3 <- read.table("caco3.txt")
caco3 <- data.frame("MYYrs"=-caco3$V1, "CaCO3"=caco3$V2)
g_ca <- caco3$MYYrs > -39 & caco3$MYYrs < -32  # Data with collapse (for plot)
p_ca <- caco3$MYYrs > -39 & caco3$MYYrs < -34  # Data used in warning signal
X <- data.frame("time"=caco3$MYYrs[p_ca], "data"=caco3$CaCO3[p_ca])
# Rather annoying to have time backwards and negative, lets reverse this.
X <- data.frame("time"=rev(X[,1] - min(X[,1])), "data"=rev(X[,2]))
## Also really annoying that there are replicate values, luckily a quick averaging call will remove them. 
require(limma)
X <-avereps(X, ID=X[,1])
pars <- c(Ro=5.0, m= -.04, theta=mean(X[,2]), sigma=sd(X[,2])*5*2)
const_pars <- c(Ro=5.0, theta=mean(X[,2]), sigma=sd(X[,2])*5*2)
const <- updateGauss(const_LTC, const_pars, X, control=list(maxit=1000))
timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))
llik_warning <- 2*(loglik(timedep)-loglik(const))

sfInit(parallel=TRUE, cpu=cpu)
sfLibrary(warningsignals)
sfExportAll()
out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
save(list=ls(), file="dakos.Rdat")
social_plot(plot(out), file="dakos.png", tag="warningsignal stochpop LTC climatedata CaCO3")





## Falsely pretend it has even time intervals
#myyrs <- caco3$MYYrs[dat]
#end=myyrs[1]; start=myyrs[length(myyrs)]
#X <- ts(caco3[dat,2],start=start,end=end, freq=length(dat)/(end-start))
#plot(caco3[dat,1], caco3[dat,2])


