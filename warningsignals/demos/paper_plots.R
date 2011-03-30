equire(warningsignals)

PNG=TRUE
PNG=FALSE

#load output of figure1.R
load("5554763401.Rdat")
source("../R/indicators.R")

indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
#indicators <- c("Variance", "Autocorrelation")


if(PNG){ png(file="figure1.png", height=length(indicators)*240, width=3*480)
} else { cairo_pdf(file="figure1.pdf",height=length(indicators)*7/2, width=3*7) }
all_indicators(	list(Deteriorating=deteriorating, Constant=constant, Empirical=deut3),	
					indicators=indicators, cex.axis=2, cex.lab=2.3)
dev.off()


## load output of figure2.R, 2000 reps
load("5554848679.Rdat")

## use just indicators 1 and 2
#deterior_taus <- deterior_taus[1:2]; constant_taus <-  constant_taus[1:2]; deut3_taus <- deut3_taus[1:2]


if(PNG) { png(file="figure2.png", width=3*480/3, height=480*length(constant_taus)/3)
} else { cairo_pdf(file="figure2.pdf", width=3*7/3, height=7*length(constant_taus)/3) }
plot.bootstrap_tau(list(Deteriorating=deterior_taus, Constant=constant_taus, Empirical=deut3_taus), 
					cex.axis=1, show_p=FALSE, ylim=c(0,2.8))
dev.off()


## load output of figure 3
#load("35555677786.Rdat")
load("35563325713.Rdat")
data_names <- c("Deteriorating", "Constant", "Empirical")

if(PNG) {  png(file="figure3.png", width=480*1.5, height=480/1.5) 
} else { cairo_pdf(file="figure3.pdf", width=7*3/3, height=3.5) }

	par(mfrow=c(1,3),  oma=c(8,8,8,4), mar=c(0,0,0,0))
	plot(deterior_mc,show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), xlim=c(0,90))
	mtext(data_names[1], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
	mtext("Probability density", WEST<-2, line=4) ## statistic name on first column

	plot(constant_mc,show_text = c("p","power"), xlab="", main="",  cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,90))
	mtext(data_names[2], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
	mtext("Likelihood Ratio", SOUTH<-1, line=4) ## x-axis label

	plot(deut3_mc,show_text = c("p","power"), xlab="", main="",  cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,90)  )
	mtext(data_names[3], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row

dev.off()




