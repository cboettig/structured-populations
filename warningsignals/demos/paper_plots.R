require(warningsignals)


#load output of figure1.R
load("5554763401.Rdat")

source("../R/indicators.R")
png(file="figure1.png", height=length(indicators)*240, width=3*480)
#cairo_pdf(file="figure1.pdf",height=length(indicators)*7/2, width=3*7)
all_indicators(	list(deteriorating=deteriorating, constant=constant, empirical=deut3),	
					indicators=indicators, cex.axis=1.8, cex.lab=1.6)
dev.off()


## load output of figure2
load("5554848679.Rdat")

png(file="figure2.png", width=3*480/3, height=480*length(constant_taus)/3)
#cairo_pdf(file="figure2.pdf", width=3*7/3, height=7*length(constant_taus)/3)
plot.bootstrap_tau(list(Deteriorating=deterior_taus, Constant=constant_taus, Empirical=deut3_taus), 
					cex.axis=1, show_p=FALSE, ylim=c(0,2.8))
dev.off()


## load output of figure 3
load("35555677786.Rdat")


png(file="figure3.png", width=480*3/2, height=480/2)
#cairo_pdf(file="figure3.pdf", width=7*3/2, height=7/2)
par(mfrow=c(1,3), mar=c(5.1,4,4.0,1),oma=c(0,0,0,0))
plot(deterior_mc,show_text = c("p","power"), xlab="Likelihood Ratio", main="Deteriorating", ylim=c(0,.7), xlim=c(0,100))
par(  mar=c(5.1,2,4.1,2) )
plot(constant_mc,show_text = c("p","power"), xlab="Likelihood Ratio", main="Constant", ylim=c(0,.7), yaxt="n", ylab="", xlim=c(0,100))
par( mar=c(5.1,0,4.1,3))
plot(deut3_mc,show_text = c("p","power"), xlab="Likelihood Ratio", main="Empirical", ylim=c(0,.7), yaxt="n", ylab="", xlim=c(0,100)  )
dev.off()


### Alternate axis setup for figure 3

png(file="figure3_alt.png", width=480*3/2, height=480/2)
#cairo_pdf(file="figure3_alt.pdf", width=7*3/2, height=7/2)
par(mfrow=c(1,3), mar=c(5.1,4,4.0,1),oma=c(0,0,0,0))
plot(deterior_mc,show_text = c("p","power"), xlab="Likelihood Ratio", main="Deteriorating", ylim=c(0,.7))
par(  mar=c(5.1,2,4.1,2) )
plot(constant_mc,show_text = c("p","power"), xlab="Likelihood Ratio", main="Constant", ylim=c(0,.7), yaxt="n", ylab="")
par( mar=c(5.1,0,4.1,3))
plot(deut3_mc,show_text = c("p","power"), xlab="Likelihood Ratio", main="Empirical", ylim=c(0,.7), yaxt="n", ylab="")
dev.off()




