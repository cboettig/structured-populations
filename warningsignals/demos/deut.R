require(socialR)
require(warningsignals)
tags<-"warningsignals stochpop climatedata deut"
cpu <- 16
nboot <- 16

load("deut_data.Rdat")


# Fit models
models <- vector("list", 4)
for(i in 1:4){
	X <- data[[i]]$X_ts
	const_pars <- c(Ro=5.0, theta=mean(X[,2]), sigma=sd(X[,2])*5*2)
	pars <- c(Ro=5.0, m= -.04, theta=mean(X[,2]), sigma=sd(X[,2])*5*2)
	const <- updateGauss(const_LTC, pars, X, control=list(maxit=1000))
	timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=1000))

	print(llik_warning <- 2*(loglik(timedep)-loglik(const)))
	models[[i]] <- list(X, const_pars, pars, const, timedep, llik_warning)
}


	sfInit(parallel=TRUE, cpu=cpu)
	sfLibrary(warningsignals)
	sfExportAll()

for(i in 1:4){
	const <- models[[i]]$const; timedep <- models[[i]]$timedep;
	out <- montecarlotest(const, timedep, cpu=cpu, nboot=nboot, GetParNames=FALSE)
	save(list=ls(), file="deut.Rdat")
	social_plot(plot(out), file="LTC_deut.png", tag=tags)
}



#plot(glaciationI)
#dim(glaciationI) ## Note this doesn't equal N=591 because they use interpolated, evenly spaced data after all




