## Various plotting funtions for simulations
plot_means <- function(beetle_data, ibm=NULL, overlay=FALSE){

	m <- max(beetle_data[,2:5])*1.2
	if(!overlay){
		lty=1
	plot(beetle_data[,1], beetle_data[,2], type = 'l', col="yellow", 
		lwd=3, ylim=c(0,m), xlab="time", ylab="mean", cex.lab=1.3, main="Beetle ELPA model" )
	} else { 
		lty=2
		lines(beetle_data[,1], beetle_data[,2], type = 'l', col="yellow", lwd=3, lty=lty)
	}
	lines(beetle_data[,1], beetle_data[,3], col="yellowgreen", lwd=3, lty=lty)	
	lines(beetle_data[,1], beetle_data[,4], col="lightgreen", lwd = 3, lty=lty)	
	lines(beetle_data[,1], beetle_data[,5], col="darkgreen", lwd = 3, lty=lty)	
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )


	if(!is.null(ibm)){
		points(beetle_data[,1], ibm$mv[[1,1]], col="yellow")	
		points(beetle_data[,1], ibm$mv[[1,2]], col="yellowgreen")	
		points(beetle_data[,1], ibm$mv[[1,3]], col="lightgreen")	
		points(beetle_data[,1], ibm$mv[[1,4]], col="darkgreen")	
	}

}


plot_stdev <- function(beetle_data, ibm=NULL, overlay=FALSE){
	v <- max(sqrt(beetle_data[,6:9]))
	if(!overlay){
		lty=1
		plot(beetle_data[,1], sqrt(beetle_data[,6]), type = 'l', col="yellow",
			lwd=3, ylim=c(0,v), xlab="time", ylab="stdev", cex=1.3 )
	} else {
		lty=2
		lines(beetle_data[,1], sqrt(beetle_data[,6]), type = 'l', col="yellow",	lwd=3, lty=lty)
	}
	lines(beetle_data[,1], sqrt(beetle_data[,7]), col="yellowgreen", lwd=3, lty=lty)	
	lines(beetle_data[,1], sqrt(beetle_data[,8]), col="lightgreen", lwd = 3, lty=lty)	
	lines(beetle_data[,1], sqrt(beetle_data[,9]), col="darkgreen", lwd = 3, lty=lty)
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )

	if(!is.null(ibm)){
		points(beetle_data[,1], sqrt(ibm$mv[[2,1]]), col="yellow")	
		points(beetle_data[,1], sqrt(ibm$mv[[2,2]]), col="yellowgreen")	
		points(beetle_data[,1], sqrt(ibm$mv[[2,3]]), col="lightgreen")	
		points(beetle_data[,1], sqrt(ibm$mv[[2,4]]), col="darkgreen")
	}
}


plot_covar <- function(beetle_data, ibm=NULL, overlay=FALSE){
	vu <- max(beetle_data[,10:15])
	vl <- min(beetle_data[,10:15])
	if(!overlay){
		lty=1
		plot(beetle_data[,1], beetle_data[,10], type = 'l', col="yellow",
			lwd=3, ylim=c(vl,vu), xlab="time", ylab="cov", cex=1.3 )
	} else {
		lty=2
		lines(beetle_data[,1], beetle_data[,10], type = 'l', col="yellow", lwd=3, lty=lty)
	}
	lines(beetle_data[,1], beetle_data[,11], col="yellowgreen", lwd=3, lty=lty)	
	lines(beetle_data[,1], beetle_data[,12], col="lightgreen", lwd = 3, lty=lty)	
	lines(beetle_data[,1], beetle_data[,13], col="darkgreen", lwd = 3, lty=lty)
	lines(beetle_data[,1], beetle_data[,14], col="lightblue", lwd = 3, lty=lty)
	lines(beetle_data[,1], beetle_data[,15], col="blue", lwd = 3, lty=lty)
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )

}



plot_stdev_div_mean <- function(beetle_data){
	eps <- .1
	v <- max(sqrt(beetle_data[,6:9])/(eps+beetle_data[,2:5]))
	plot(beetle_data[,1], sqrt(beetle_data[,6])/beetle_data[,2], type = 'l', col="yellow",
		lwd=3, ylim=c(0,v), xlab="time", ylab="stdev/mean", cex=1.3 )
	lines(beetle_data[,1], sqrt(beetle_data[,7])/beetle_data[,3], col="yellowgreen", lwd=3)	
	lines(beetle_data[,1], sqrt(beetle_data[,8])/beetle_data[,4], col="lightgreen", lwd = 3)	
	lines(beetle_data[,1], sqrt(beetle_data[,9])/beetle_data[,5], col="darkgreen", lwd = 3)
	legend("right", c("egg", "larva", "pupa", "adult"), 
		lty=1, col=c("yellow", "yellowgreen", "lightgreen", "darkgreen") )
}




