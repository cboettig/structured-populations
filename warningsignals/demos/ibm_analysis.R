#ibm_analysis.R

load("ibm_sims.Rdat")
indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
nboot <- 80
cpu <- 16
ibm_tags="warningsignals stochpop ibm", 

social_plot(
	all_indicators(	list(deteriorating=ibm_critical, 
						 constant=ibm_stable),	
					indicators=indicators, method="kendall"
	),
	tags=ibm_tags, 
	height=length(indicators)*240, width=2*480
)


## fit models
deterior_m<-fit_models(ibm_critical, "LSN")
constant_m<-fit_models(ibm_stable, "LSN")
deterior_taus <- bootstrap_tau(deterior_m$X,
							   deterior_m$const, deterior_m$timedep,
							   indicators = indicators,
							   method="kendall",
							   nboot=nboot, cpu=cpu)
constant_taus <- bootstrap_tau(constant_m$X, constant_m$const,
							   constant_m$timedep,
							   indicators = indicators,
							   method="kendall",
							   nboot=nboot, cpu=cpu)
social_plot(
	plot.bootstrap_tau(list(deteriorating=deterior_taus,
		constant=constant_taus),
		show_p=TRUE, cex.axis=2),
	tags = ibm_tags, height=480*2, width=480*2
	)

deterior_mc <- 
		montecarlotest(deterior_m$const, deterior_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)

constant_mc <- 
		montecarlotest(constant_m$const, constant_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)



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

modelchoice_mc <- 
		montecarlotest(deterior_LTC$timedep, deterior_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)
social_plot(plot(modelchoice_mc), tags=ibm_tags)

