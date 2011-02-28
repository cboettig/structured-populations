# data from: 
deut <- read.table("deutnat.txt")

# divisions from: Dakos et. al.
g1 <- which(-deut$V2 > -58800 & -deut$V2 < -12000)
p1 <- which(-deut$V2 > -58800 & -deut$V2 < -17000)
glaciationI   <- data.frame("time"=-deut$V2[p1], "data"=deut$V3[p1])
g2 <- which(-deut$V2 > -151000 & -deut$V2 < -128000)
p2 <- which(-deut$V2 > -151000 & -deut$V2 < -135000)
glaciationII   <- data.frame("time"=-deut$V2[p2], "data"=deut$V3[p2])
g3 <- which(-deut$V2 > -270000 & -deut$V2 < -238000)
p3 <- which(-deut$V2 > -270000 & -deut$V2 < -242000)
glaciationIII   <- data.frame("time"=-deut$V2[p3], "data"=deut$V3[p3])
g4 <- which(-deut$V2 > -385300 & -deut$V2 < -324600)
p4 <- which(-deut$V2 > -385300 & -deut$V2 < -334100)
glaciationIV   <- data.frame("time"=-deut$V2[p4], "data"=deut$V3[p4])


X <- glaciationI
pars <- c(Ro=5.0, m= -.04, theta=mean(X[,2]), sigma=sd(X[,2]))
const_pars <- c(Ro=5.0, theta=mean(X[,2]), sigma=sd(X[,2]))
timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))
const <- updateGauss(const_LTC, pars, X, control=list(maxit=1000))
llik_warning <- 2*(loglik(timedep)-loglik(const))


plot(glaciationI)
dim(glaciationI) ## Note this doesn't equal N=591 because they use interpolated, evenly spaced data after all



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



