source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
##### WHY ARE THESE NECESSARY???
load("5555038554.Rdat")
source("../R/indicators.R")


indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
## FIGURE 1
social_plot(
	all_indicators(list(CaCO3=X, 
				        GlaciationI=data[[1]]$X_ts, 
				        GlaciationII=data[[2]]$X_ts ),	
				   indicators=indicators, method="kendall"),
	tags="warningsignals stochpop", 
	height=length(indicators)*240, width=3*480
)

cairo_pdf(file="figure1_appendix.pdf",height=7*length(indicators)/2, width=3*7)
all_indicators(list(CaCO3=X, 
				    GlaciationI=data[[1]]$X_ts, 
				    GlaciationII=data[[2]]$X_ts ),	
				indicators=indicators, method="kendall"
                cex.axis=2, cex.lab=2.3)
dev.off()




## load output of figure2.R, for the CaCO3 data
load("5554848679.Rdat")
## Load model fits and tau bootstraps.  This file produced by deut_examples.Rdat
load("5562383846.Rdat")
social_plot(
	plot.bootstrap_tau(
		list(CaCO3=CaCO3_taus, GlaciationI=taus[[1]], GlaciationII=taus[[2]], 
		cex.axis=1, ylim = c(0,2.8)), tags = "warningsingals stochpop", 
			height=480*2, width=480*3
	)


## load output of figure3.R, for the CaCO3 data
load("35555677786.Rdat")
## load the likelihood bootstraps, produced by deut_examples.Rdat -- should rerun with more replicates
load("5562961240.Rdat")

deut_labels <- c("GlaciationI", "GlaciationII")

plt <- function(){
	par(mfrow=c(1,3),  oma=c(4,4,4,4), mar=c(0,0,0,0))
	plot(mc[[i]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4))
	mtext("CaCO3", NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
	for(i in 1:2){
		plot(mc[[i]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4))
   		mtext(deut_labels[i], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
 }
		mtext("Probability density", WEST<-2, outer=TRUE, line=2) ## statistic name on first column
		mtext("Likelihood Ratio", SOUTH<-1, outer=TRUE, line=2) ## x-axis label
}

social_plot(plt(), 	tags = "warningsingals stochpop", 
			height=480, width=480*3)


