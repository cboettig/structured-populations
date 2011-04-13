# montecarlotest.R
## Dependencies: snowfall (parallel computing)


## requries loglik and getParameters to be defined functions. Here we just set up there generics.   

# loglik isn't a standard S3 generic, so we add it:
loglik <- function(x, y, ...) UseMethod("loglik")

# getParameters isn't a standard S3 generic, so we add it:
getParameters <- function(x, y, ...) UseMethod("getParameters")


### Threshold should really be handled by the plotting function...
montecarlotest <- function(null, test, nboot = 100, cpu = 2, 
	threshold = .95, GetParNames=TRUE, method= c("Nelder-Mead", 
	"BFGS", "CG", "L-BFGS-B", "SANN"), ...){
	## are we in parallel?
	if(cpu>1 & !sfIsRunning()){ 	
		sfInit(parallel=TRUE, cpu=cpu) 
		sfLibrary(warningsignals) ## needs to export the library with all function defs!
		sfExportAll()
	} else if(cpu<2 & !sfIsRunning()){  sfInit()
	} else { }

	## generate each distribution
	null_sim <- sfSapply(1:nboot, function(i){
		data <- simulate(null)
		if (is(data,"list")) data <- data$rep.1 # handle data that has replicates
		null <- update(null, data, method=method, ...)
		test <- update(test, data, method=method, ...)
		lr <- -2*(loglik(null) - loglik(test)) 
		list(lr, getParameters(null), getParameters(test))
	})
	test_sim <- sfSapply(1:nboot, function(i){
		data <- simulate(test)
		if (is(data,"list")) data <- data$rep.1 # handle data that has replicates
		null <- update(null, data, method=method, ...)
		test <- update(test, data, method=method, ...)
		lr <- -2*(loglik(null) - loglik(test))
		list(lr, getParameters(null), getParameters(test))
	})
	## grab the distribution of likelihood ratios
	null_dist <- unlist(null_sim[1,]) 
	test_dist <- unlist(test_sim[1,]) 

	
	## grab distribution of parameters observed for null model under null model sims
	null_bootstrap_pars <- matrix( unlist(null_sim[2,]), ncol=nboot)
	if(GetParNames) rownames(null_bootstrap_pars) <- names(getParameters(null))

	## grab distribution of parameters observed for test model under test model sims
	test_bootstrap_pars <- matrix( unlist(test_sim[3,]), ncol=nboot)
	if(GetParNames) rownames(test_bootstrap_pars) <- names(getParameters(test))

	### Grab the cross-case bootstraps too:
	## grab distribution of parameters observed for test model under null model sims
	null_sim_test_pars <- matrix( unlist(null_sim[3,]), ncol=nboot)
	if(GetParNames) rownames(null_sim_test_pars) <- names(getParameters(test))

	## grab distribution of parameters observed for null model under test model sims
	test_sim_null_pars <- matrix( unlist(test_sim[2,]), ncol=nboot)
	if(GetParNames) rownames(test_sim_null_pars) <- names(getParameters(null))


	## are the parameter distribution estimates of one model under the other' simulations interesting?

	## Power calculation
	threshold_tail <- sort(null_dist)[ round(threshold*nboot) ]
	power <- sum(test_dist > threshold_tail)/nboot

	lr <- -2*(loglik(null) - loglik(test))
	p <- 1-sum(null_dist < lr)/nboot
	print(paste("power is ", power, ", p = ", p))

	# Check the reverse test (equivalent to swapping null and test)
	reverse_p <- sum(test_dist < lr)/nboot
	reverse_threshold_tail <- sort(-test_dist)[ round(threshold*nboot) ]
	reverse_power <- sum(-null_dist > reverse_threshold_tail)/nboot
	print(paste("reverse test power ", reverse_power, ", reverse test p = ", reverse_p))


	## format the output
	output <- list(power=power, null_dist=null_dist, test_dist=test_dist, nboot=nboot, 
					threshold=threshold, null=null, test=test, p = p, lr = lr, 
					reverse_p = reverse_p, reverse_power = reverse_power, 
					null_par_dist = null_bootstrap_pars, test_par_dist = test_bootstrap_pars,
					null_sim_test_pars = null_sim_test_pars, test_sim_null_pars=test_sim_null_pars)
	class(output) <- "pow"
	output
}

KL_divergence <- function(null_dist, test_dist){
	nd <- density(null_dist)
	td <- density(test_dist)
	P <- nd$y[nd$y > 0 & td$y > 0]
	Q <- td$y[nd$y > 0 & td$y > 0]
	KL_div <- P %*% log(P/Q)
	KL_div	
}

overlap <- function(pow, bw="nrd0"){
	nd <- density(pow$null_dist, bw=bw)
	td <- density(pow$test_dist, bw=bw)
	td$y %*% nd$y
}


## can now specify the info criterion
## plotting function
plot.pow <- function(pow, main="", legend=FALSE, type="density", test_dist=TRUE, shade_power=FALSE, shade_p=FALSE, show_aic=FALSE, show_data=TRUE, shade=TRUE, shade_aic=FALSE, print_text=TRUE, show_text = c("p", "power", "reverse_p", "aic"), xlim=NULL, ylim=NULL, null_dist=TRUE, bw = "nrd0", info_criterion=c("aic", "bic", "aicc", "threshold"), ...){


	## DOF calculation ##!! Should be made into a generic!!
	dof <- function(object){
		if(is(object, "fitContinuous")) dof<- object[[1]]$k
		else if(is(object, "hansentree")){
			dof <- length(object@sqrt.alpha)+length(object@sigma)+sum(sapply(object@theta,length))
		} else if(is(object, "browntree")) { 
			dof <- length(object@sigma)+sum(sapply(object@theta,length))
		} else {
			dof <- object$k
			if(is.null(object$k)) print(paste("cannot determine degrees of freedom, please input for model"))
		}
		dof
	}

	## Calculate the densities
	nd <- density(pow$null_dist, bw=bw)
	td <- density(pow$test_dist, bw=bw)


	## Calculate Axis Limits
	if(is.null(xlim)) xlim <- c( min(pow$null_dist, pow$test_dist), max(pow$null_dist, pow$test_dist) )
	if(is.null(ylim)) ylim <- c(min(td$y, nd$y), max(td$y,nd$y))


	## Select the information criterion
	info_criterion = match.arg(info_criterion)
	if(info_criterion=="aic"){
		aic_line <-  2*(dof(pow$test) - dof(pow$null)) 
	} else if(info_criterion=="threshold") {
		threshold_tail <- sort(pow$null_dist)[ round(pow$threshold*pow$nboot) ]
		aic_line <- threshold_tail #nd$x[tmp]
		print(paste("threshold", aic_line))
	}

	## Density plots
	if(type != "hist"){
        if(shade_aic==FALSE){ 
            plottype="n"
        } else { 
            plottype ="s"
        }
		## Plot the null distribution with appropriate shading
		if(null_dist){ 
			plot(nd, xlim=xlim, ylim=ylim, main=main, type=plottype, col=rgb(0,0,1,1), ...) 
			if(shade_p){
				shade_p <- which(nd$x > pow$lr)
				polygon(c(pow$lr,nd$x[shade_p]), c(0,nd$y[shade_p]), col=rgb(0,0,1,.3), border=rgb(0,0,1,.5))
			} else if(shade){
				polygon(nd$x, nd$y, col=rgb(0,0,1,.3), border=rgb(0,0,1,.5))
			}
		} else { 
			plot(0,0, xlim=xlim, ylim = ylim, main=main, xlab=" Likelihood Ratio", ...)  ## blank plot
		}

		## Plot the test distribution with appropriate shading
		if(test_dist){
			if(!null_dist) plot(td, xlim=xlim, main=main, type=plottype, col=rgb(1,0,0,1), ...)  ## just plot test dist 
			else lines(td, type=plottype, col=rgb(1,0,0,1))
			threshold_tail <- sort(pow$null_dist)[ round(pow$threshold*pow$nboot) ]
			if(shade_power){
				shade_power <- which(td$x > threshold_tail)
				polygon(c(threshold_tail, td$x[shade_power]), c(0,td$y[shade_power]), col=rgb(1,0,0,.3), border=rgb(1,0,0,.5))
			} else if(shade){
				polygon(td$x, td$y, col=rgb(1,0,0,.3), border=rgb(1,0,0,.5))
			}
		}
		## AIC shading 
		if(shade_aic){
			shade_p <- which(nd$x > aic_line)
			polygon(c(aic_line,nd$x[shade_p]), c(0,nd$y[shade_p]), col=rgb(1,.5,0,.5), border=rgb(0,0,1,0))
			if(test_dist){ # If test_dist on, shade the errors under the test too
				shade_p <- which(td$x < aic_line)
				polygon(c(aic_line,td$x[shade_p]), c(0,td$y[shade_p]), col=rgb(1,1,0,.5), border=rgb(1,0,0,0))
			}

		}

	## Plot histograms instead of density plots
	} else {
	hist(pow$null_dist, xlim=xlim, lwd=3, col=rgb(0,0,1,.5), border="white", main=main, ...)
		if(test_dist){
			hist(pow$test_dist, add=T, lwd=0, col=rgb(1,0,0,.5), border="white", main=main, ...)
		}
	}

	##  Add lines to plots 
	if(show_data){ 
		#abline(v=pow$lr, lwd=3, col="darkred", lty=2 )
		points(pow$lr,yshift(1), cex=1.5, col="black", pch=25, fg="black", bg="black")
	}
	if(show_aic){
    #    abline(v=aic_line, lwd=3, col="darkgray", lty=3) 
    	points(aic_line,yshift(1), cex=1.5, col="red", pch=25, fg="red", bg="red")
}

	## Calculate statistics
	p <- pow$p
	aic_wrong <- sum(pow$null_dist > aic_line)/pow$nboot
	aic_power <- sum(pow$test_dist > aic_line)/pow$nboot





	## Print statistics to console, can display on plot
#	if(aic_line < pow$lr){ print("AIC prefers test model")} else { print("AIC prefers null model") }	
	print(paste("p = ",  p, ", power = ", pow$power))
	if(print_text){
		if (! is.na(match("p", show_text)) )
            text(xshift(105), yshift(85), paste("Type I = ", round(p,3), 
            "\n Type II = ", round(1-pow$power,3)),pos=2, cex=par()$cex.lab/.8)
#		if (! is.na(match("power", show_text)) )
 #           text(xshift(105), yshift(85), paste("Type II = ", round(1-pow$power,3)), pos=2)
		if (! is.na(match("aic", show_text)) )
            text(xshift(105), yshift(75), paste("AIC false alarm rate", aic_wrong*100, "%"), pos=2)
		if (! is.na(match("aic", show_text)) )
            text(xshift(105), yshift(70), paste("AIC missed detection rate", (1-aic_power)*100, "%"), pos=2)
    	if (! is.na(match("reverse_p", show_text)) )
            text(xshift(105), yshift(65), paste("p in reverse test ", pow$reverse_p), pos=2)
	}
	## add legend
	if(legend){
        if (shade==TRUE){
		    legend("topright", c("sim under Null", "sim under Test", "observed"), 
                   pch=c(15,15,25), fg=c("white", "white", "black"), col=c(rgb(0,0,1,.5), rgb(1,0,0,.5), "black"))
        }
        else if (shade_aic==TRUE & test_dist==TRUE){
   		    legend("topright", c( paste("False Alarms (", round(aic_wrong*100,3), "%)", sep=""),
                                  paste("Missed Events (", round((1-aic_power)*100,3), "%)", sep="")), 
                                  pch=c(15,15), col=c(rgb(1,.5,0,.5), rgb(1,1,0,.5)))
        }
        else if (shade_aic==TRUE & test_dist==FALSE){
   		    legend("topright", paste("False Alarms (", round(aic_wrong*100,3), "%)", sep=""),
                                  pch=15, col=rgb(1,.5,0,.5))
        }
    }
}

