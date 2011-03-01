require(socialR)
require(warningsignals)
tags<-"warningsignal stochpop LTC climatedata deut"
cpu <- 16
nboot <- 160


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

glaciation <- list(glaciationI, glaciationII, glaciationIII, glaciationIV)
for(i in 1:4){
	X <- glaciation[[i]]
	require(limma)
	X <-avereps(X, ID=X[,1])
	pars <- c(Ro=5.0, m= -.04, theta=mean(X[,2]), sigma=sd(X[,2])*5*2)
	const_pars <- c(Ro=5.0, theta=mean(X[,2]), sigma=sd(X[,2])*5*2)
	timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))
	const <- updateGauss(const_LTC, pars, X, control=list(maxit=1000))

	print(llik_warning <- 2*(loglik(timedep)-loglik(const)))
	sfInit(parallel=TRUE, cpu=cpu)
	sfLibrary(warningsignals)
	sfExportAll()

	out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
	save(list=ls(), file="deut.Rdat")
	social_plot(plot(out), file="LTC_deut.png", tag=tags)
}



#plot(glaciationI)
#dim(glaciationI) ## Note this doesn't equal N=591 because they use interpolated, evenly spaced data after all




