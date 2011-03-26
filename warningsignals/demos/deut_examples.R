## FIGURE 2
require(warningsignals)
require(socialR)
gitcommit()

## load datafiles
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts
load("5555038554.Rdat") # an LSN sim, stable and collapsing
#load("5550863240.Rdat") # an LTC sim, just collapsing
#simLTC <- X

indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
nboot <- 80
cpu <- 16


m <- lapply(1:4, function(i){
	X <- data[[i]]$X_ts
	fit_models(X, "LSN")
})

taus <- lapply(1:4, function(i){
	bootstrap_tau(m[[1]]$X, m[[1]]$const, m[[1]]$timedep, indicators=indicators, nboot=nboot, cpu=cpu)
})

social_plot(
	plot.bootstrap_tau(list(I=taus[[1]], II=taus[[2]], III=taus[[3]], IV=taus[[4]]), cex.axis=2),
			tags = "warningsingals stochpop", 
			height=480*2, width=480*4
	)

mc <- lapply(1:4, function(i){
	montecarlotest(m[[1]]$const, m[[1]]$timedep, nboot=nboot, cpu=cpu, GetParNames=FALSE)
})

save(list=ls(), file="deut_examples.Rdat")

plt <- function(){
	par(mfrow=c(1,3),  oma=c(8,8,8,4), mar=c(0,0,0,0))
	for(i in 1:4){
		plot(mc[[i]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), xlim=c(0,90))
		mtext("Probability density", WEST<-2, outer=TRUE) ## statistic name on first column
		mtext(paste(i), NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
		mtext("Likelihood Ratio", SOUTH<-1, outer=TRUE) ## x-axis label
	}
}

social_plot(plt(), 	tags = "warningsingals stochpop", 
			height=480*1.5, width=480*4)


