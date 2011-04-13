# bootstrap_indicators.R

## A few wrapper functions to make it easy to bootstrap a secified set of indicator statistics

fit_models <- function(X, model=c("LTC", "LSN"), integrateOU=FALSE,  
					   optim_method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
					   ...){
	optim_method <- match.arg(optim_method)
	lower <- -Inf
	upper <- Inf
	if( optim_method=="L-BFGS-B"){
		lower = c(0, -Inf, -Inf, 0)
		upper = c(Inf, 0, Inf, Inf)
	}
	const_pars <- c(Ro=as.numeric(1/max(time(X))), theta=as.numeric(mean(X)), sigma=as.numeric(sd(X)))
# Fit a linearized transcritical bifurcation model
	if (integrateOU){
		if (model=="LSN"){
			const <- updateGauss(const_LSN, const_pars, X, method=optim_method, 
								 control=list(maxit=2000), ...)
		} 
		else if (model=="LTC"){
			const <- updateGauss(const_LTC, const_pars, X, method=optim_method, 
								 control=list(maxit=2000), ...)
		}
	} 
	else {	
	const <- updateGauss(constOU, const_pars, X, method=optim_method, 
						 control=list(maxit=2000), ...)
	}

	if (model=="LTC"){
    	pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), 
    			  sigma=as.numeric(const$pars["sigma"]))
		timedep <- updateGauss(timedep_LTC, pars, X, method=optim_method, 
							   control=list(maxit=2000), upper=upper, lower=lower, ...)
	} 
	else if (model=="LSN"){

        # guess LSN parameters from OU parameterization
        guess_Ro <- as.numeric(const$pars['Ro']^2)
        guess_theta <- as.numeric(const$pars['theta']+const$pars['Ro'] )
        guess_sigma<- as.numeric(const$pars['sigma']/sqrt(2*const$pars['Ro']+const$pars['theta']))
    	pars <- c(Ro=guess_Ro, m=0, theta=guess_theta, sigma=guess_sigma)

		timedep <- updateGauss(timedep_LSN, pars, X, method=optim_method, 
							   control=list(maxit=2000), upper=upper, lower=lower, ...)
	}
	list(X=X, const=const, timedep=timedep, pars=pars, const_pars=const_pars,
		 model=model)
}


bootstrap_tau <- function(X, const, timedep, 
						  indicators = c("Variance", "Autocorrelation", "Skew", "Kurtosis", "CV"),
						  nboot=160, cpu=16, windowsize=round(length(X)/2), 
						  method=c("pearson", "kendall", "spearman")){
# Tau approach comparison
	taus <- lapply(indicators, 
				   function(stat){
						tau_dist_montecarlo(X, const, timedep, signal=stat, 
						                    nboot=nboot, cpu=cpu, method=method) 
						})
	class(taus) <- "bootstrap_tau"
	taus
}


plot.bootstrap_tau <- function(taus, show_p = FALSE, show_error=TRUE, ...){
## If elements of "taus" are of class tau_dist_montecarlo, assumes we have 
## bootstraps for only a single dataset, and plot a single column.   
## If elements are also "bootstrap_taus" then we have a data.frame 
## with bootstraps from multiple datasets, and we plot a matrix with column 
## for each dataset and row for each indicator
	if(is(taus[[1]], "tau_dist_montecarlo")){
	## treat single dataset as data.frame notation as well
		taus <- list(taus)
	}
	data_names <- names(taus)

	## dimensions
	n <- length(taus) ## number of datasets
	m <- length(taus[[1]]) ## number of indicators

	## set up m x n plot-matrix, no margins on subplots, add outer margins
	par(mfrow=c(m,n), oma=c(3,4,2,.2), mar=c(0,0,0,0), ...)

	for(j in 1:m){
		for(i in 1:n){
			if(j == m){ xaxt <- "s"
			} else { xaxt <- "n" }
			if(i > 1){ yaxt <- "n" 
			} else { yaxt <- "s" }

			plot(taus[[i]][[j]], show_p = show_p, show_error=show_error, 
                 xaxt = xaxt, yaxt=yaxt, ...)
			if(j==1) mtext(data_names[i], NORTH<-3, 
                           cex=par()$cex.lab, line=1) ## data labels on top
			if(i==1){
				mtext(taus[[i]][[j]]$signal, WEST<-2, line=3, cex=par()$cex.lab) ## statistic name on first column
				mtext(expression(paste("Prob Density of ", tau)),
						WEST<-2, line=2, cex=.8*par()$cex.lab) ## statistic name 
			}	
			if(j==m & i==2) mtext(expression(paste("Correlation coefficient, ", tau)), SOUTH<-1, line=2, cex=par()$cex.lab) ## x-axis label
		}
	}
}

