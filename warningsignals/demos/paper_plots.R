require(warningsignals)

JPEG=TRUE

#load output of figure1.R
load("5554763401.Rdat")
source("load_deut.R")
source("../R/indicators.R")
source("../R/bootstrap_indicators.R")

#indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
indicators <- c("Variance", "Autocorrelation")


if(JPEG){ jpeg(file="Boettiger_fig1.jpg", height=length(indicators)*183/4*(3/5), width=183, units="mm", quality=100, res=150)
} else { cairo_pdf(file="figure1.pdf",height=length(indicators)*7/2, width=3*7) }
all_indicators(	list(Deteriorating=deteriorating, Constant=constant, 
                GlaciationI=data[[1]]$X_ts, GlaciationII=data[[2]]$X_ts, GlaciationIII=data[[3]]$X_ts),	
				indicators=indicators, cex.axis=.5, cex.lab=.6, lwd=.5, yaxs="i")
dev.off()


## load output of figure2.R, 2000 reps
load("5554848679.Rdat")
## Load model fits and tau bootstraps.  This file produced by deut_examples.Rdat with 2000 replicates
load("5562383846.Rdat")  #kendall, 2000

## use just indicators 1 and 2: as in deterior_taus[1:2] (var and autocorr)

if (JPEG){ jpeg(file="Boettiger_fig2.jpg", height=2*37, width=183, units="mm", quality=100, res=150)
} else { cairo_pdf(file="figure2.pdf", width=3*7/3, height=7*length(constant_taus)/3) }
plot.bootstrap_tau(list(Deteriorating=deterior_taus[1:2], 
                   Constant=constant_taus[1:2], GlaciationI=taus[[1]][1:2], 
                   GlaciationII=taus[[2]][1:2], GlaciationIII=taus[[3]][1:2]), 
				   cex.axis=.6, cex.lab=.8, show_p=FALSE, ylim=c(0,2.8), yaxp = c(0, 3, 3), xaxp=c(-1,1,5))
dev.off()


## load output of figure 3
load("35563325713.Rdat")
## 2000 replicates from deut_likelihood.R
load("5592395409.Rdat")


data_names <- c("Deteriorating", "Constant", "Empirical")

deut_labels <- c("GlaciationI", "GlaciationII", "GlaciationIII")

if(JPEG) { jpeg(file="Boettiger_fig3.jpg", height=37*1.2, width=183, units="mm", quality=100, res=150)
} else { cairo_pdf(file="figure3.pdf", width=7*3/3, height=3.5) }

	par(mfrow=c(1,5),  oma=c(3,3,3,1), mar=c(0,0,0,0), cex.lab=.8, cex.axis=.8)
	plot(deterior_mc,show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), xlim=c(0,40))
	mtext(data_names[1], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row
	mtext("Probability density", WEST<-2, line=2, cex=par()$cex.lab) ## statistic name on first column

	plot(constant_mc,show_text = c("p","power"), xlab="", main="",  cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,40))
	mtext(data_names[2], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row
#	mtext("Likelihood Ratio", SOUTH<-1, line=2, cex=par()$cex.lab) ## x-axis label

#	plot(deut3_mc,show_text = c("p","power"), xlab="", main="",  cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=c(0,90)  )
#	mtext(data_names[3], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row
    xlim <- c(0,600)
	for (i in 1:3){
        if(i>1) xlim <- c(0,100)
		plot(mc[[i]],show_text = c("p","power"), xlab="", main="", cex.lab=1, ylim=c(0,.4), yaxt="n", ylab="", xlim=xlim)
   		mtext(deut_labels[i], NORTH<-3, cex=par()$cex.lab, line=1) ## data labels on top row
        if (i==1)
            mtext("Likelihood Ratio", SOUTH<-1, line=2, cex=par()$cex.lab) ## x-axis label
}
	
dev.off()




