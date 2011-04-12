require(warningsignals)


#load output of figure1.R
load("5554763401.Rdat")  ## Simulated data from LSN with and without an impending collapse.  
deterior <- deteriorating  ## Watch out for name collisions with all these data sets
source("load_deut.R")    ## Deuterium data runs
## should load the CaCO3 data and assign it to a name if you want it

## Drake Data example, for all figure types
load("reanalyze.Rdat")
drake <- out

source("../R/indicators.R")
source("../R/bootstrap_indicators.R")


jpeg(file="Boettiger_fig1.jpg", height=length(indicators)*183/4*(3/5), width=183, units="mm", quality=100, res=150)
all_indicators(	list(Constant=constant, Deteriorating=deterior, 
                Glaciation=data[[3]]$X_ts, Algae=drake$data[["H6"]]),	
				indicators=indicators, cex.axis=.5, cex.lab=.6, lwd=.5, yaxs="i")
dev.off()

### APPENDIX 
appendix_indicators <- c("Skew", "Kurtosis")
jpeg(file="Boettiger_fig1_appendix.jpg", height=length(appendix_indicators)*183/4*(3/5), width=183, units="mm", quality=100, res=150)
all_indicators(	list(Constant=constant, Deteriorating=deterior, 
                Glaciation=data[[3]]$X_ts, Algae=drake$data[["H6"]]),	
				indicators=appendix_indicators, cex.axis=.5, cex.lab=.6, lwd=.5, yaxs="i")
dev.off()


## load output of figure2.R, 2000 reps
# contains the simulated data from LSN (deterior and constant), deut3 and CaCO3 data
load("5554848679.Rdat") 
## Load model fits and tau bootstraps.  This file produced by deut_examples.Rdat with 2000 replicates
## Contains deuterium data from Glaciation I, II, III using Kendall's test
load("5562383846.Rdat")  #kendall, 2000

## use just indicators 1 and 2: as in deterior_taus[1:2] (var and autocorr)
jpeg(file="Boettiger_fig2.jpg", height=2*37, width=183, units="mm", quality=100, res=150)
plot.bootstrap_tau(list(Constant=constant_taus[1:2],
                   Deteriorating=deterior_taus[1:2], 
                   Glaciation=taus[[3]][1:2], 
                   Algae=drake$taus[["H6"]][1:2]),
				   cex.axis=.6, cex.lab=.8, show_p=FALSE, ylim=c(0,2.8), yaxp = c(0, 3, 3), xaxp=c(-1,1,5))
dev.off()

########## APPENDIX PLOT
## use just indicators 3 and 4: skew & Kurtosis
jpeg(file="Boettiger_fig2_appendix.jpg", height=2*37, width=183, units="mm", quality=100, res=150)
plot.bootstrap_tau(list(Constant=constant_taus[3:4],
                   Deteriorating=deterior_taus[3:4], 
                   Glaciation=taus[[3]][3:4], 
                   Algae=drake$taus[["H6"]][3:4]),
				   cex.axis=.6, cex.lab=.8, show_p=FALSE, ylim=c(0,2.8), yaxp = c(0, 3, 3), xaxp=c(-1,1,5))
dev.off()


## load output of Figure 3: deterior (LSN sim), constant (OU sim), deut3, CaCO3 data
load("35563325713.Rdat")
## 2000 replicates from deut_likelihood.R, data on: GlaciationI, II, III
load("5592395409.Rdat")


data_names <- c( "Constant", "Deteriorating", "Glaciation", "Algae")


jpeg(file="Boettiger_fig3.jpg", height=37*1.2, width=183, units="mm", quality=100, res=150)
	par(mfrow=c(1,length(data_names)),  oma=c(3,3,2,.2), mar=c(0,0,0,0), cex.lab=.8, cex.axis=.8)
	plot(constant_mc,show_text = c("p","power"), xlab="", main="",  cex.lab=1, ylim=c(0,.4), xlim=c(0,40))
	mtext(data_names[1], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row
	mtext("Probability density", WEST<-2, line=2, cex=par()$cex.lab) ## statistic name on first column

    plot(deterior_mc,show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,40))
	mtext(data_names[2], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row

    mtext("Likelihood Ratio", SOUTH<-1, line=2, cex=par()$cex.lab) ## x-axis label

	plot(mc[[3]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,90))
   	mtext(data_names[3], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row

   	plot(drake$mc[["H6"]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,40))
   	mtext(data_names[4], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row

   	
dev.off()


require(socialR)
#social_report(file="Boettiger_fig1.jpg", tags="figure1 stochpop warningsignals publish", public=0)
#social_report(file="Boettiger_fig2.jpg", tags="figure2 stochpop warningsignals publish", public=0)
#social_report(file="Boettiger_fig3.jpg", tags="figure3 stochpop warningsignals publish", public=0)
social_report(file="Boettiger_fig1_appendix.jpg", tags="figure1 appendix stochpop warningsignals publish", public=0)
social_report(file="Boettiger_fig2_appendix.jpg", tags="figure2 appendix stochpop warningsignals publish", public=0)


