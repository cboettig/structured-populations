load("5571632581.Rdat")
require(warningsignals)
require(socialR)


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


