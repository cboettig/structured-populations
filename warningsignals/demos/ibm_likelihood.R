load("5572028881.Rdat")
require(socialR)
require(warningsignals)
indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
nboot <- 2000
cpu <-16 
ibm_tags="warningsignals stochpop ibm" 

deterior_mc <- 
		montecarlotest(deterior_m$const, deterior_m$timedep, 
		cpu=cpu, nboot=nboot)

save(list=ls(), file="ibm_likelihood.Rdat")
constant_mc <- 
		montecarlotest(constant_m$const, constant_m$timedep, 
		cpu=cpu, nboot=nboot)


save(list=ls(), file="ibm_likelihood.Rdat")

plt <- function()
	par(mfrow=c(1,2),  oma=c(6,6,6,4), mar=c(0,0,0,0))
	plot(deterior_mc,show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4))
	mtext(data_names[1], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
	mtext("Probability density", WEST<-2, line=4) ## statistic name on first column

	plot(constant_mc,show_text = c("p","power"), xlab="", main="",  cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="")
	mtext(data_names[2], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
	mtext("Likelihood Ratio", SOUTH<-1, line=4) ## x-axis label
	mtext(data_names[3], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row

social_plot(plt(), tags=ibm_tags)



deterior_LTC <-fit_models(ibm_critical, "LTC")
save(list=ls(), file="ibm_likelihood.Rdat")

modelchoice_mc <- 
		montecarlotest(deterior_LTC$timedep, deterior_m$timedep, 
		cpu=cpu, nboot=nboot)
social_plot(plot(modelchoice_mc), tags=ibm_tags)

