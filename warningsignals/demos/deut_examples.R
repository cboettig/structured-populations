## FIGURE 2
require(warningsignals)
require(socialR)
gitcommit()
## load datafiles
source("load_deut.R")

indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
nboot <- 2000
cpu <- 16


m <- lapply(1:3, function(i){
	X <- data[[i]]$X_ts
	fit_models(X, "LSN")
})

save(list=ls(), file="deut_examples.Rdat")

taus <- lapply(1:3, function(i){
	bootstrap_tau(m[[i]]$X, m[[i]]$const, m[[i]]$timedep, indicators=indicators, nboot=nboot, cpu=cpu)
})

social_plot(
	plot.bootstrap_tau(
		list(I=taus[[1]], II=taus[[2]], III=taus[[3]]), 
		cex.axis=2, ylim = c(0,2)), tags = "warningsingals stochpop", 
			height=480*2, width=480*3
	)

mc <- lapply(1:3, function(i){
	montecarlotest(m[[i]]$const, m[[i]]$timedep, nboot=nboot, cpu=cpu, GetParNames=FALSE)
})

save(list=ls(), file="deut_examples.Rdat")

deut_labels <- c("GlaciationI", "GlaciationII", "GlaciationIII")


plt <- function(){
	par(mfrow=c(1,3),  oma=c(4,4,4,4), mar=c(0,0,0,0))
	for(i in 1:3){
		plot(mc[[i]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4))
   		mtext(deut_labels[i], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
 }
		mtext("Probability density", WEST<-2, outer=TRUE, line=2) ## statistic name on first column
		mtext("Likelihood Ratio", SOUTH<-1, outer=TRUE, line=2) ## x-axis label
}

social_plot(plt(), 	tags = "warningsingals stochpop", 
			height=480/1.5, width=480*1.5)


