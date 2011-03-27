# bootstrap_indicators.R

## A few wrapper functions to make it easy to bootstrap a secified set of indicator statistics

fit_models <- function(X, method=c("LTC", "LSN")){
	const_pars <- c(Ro=1/max(time(X)), theta=mean(X), sigma=sd(X))
## Fit a linearized transcritical bifurcation model
	const <- updateGauss(constOU, const_pars, X, control=list(maxit=2000))
	pars <- c(Ro=as.numeric(const$pars["Ro"]), m=0, theta=mean(X), sigma=as.numeric(const$pars["sigma"]))
	if(method=="LTC"){
		timedep <- updateGauss(timedep_LTC, pars, X, control=list(maxit=2000))
	} else if(method=="LSN"){
		timedep <- updateGauss(timedep_LSN, pars, X, control=list(maxit=2000))
	}
	list(X=X, const=const, timedep=timedep, pars=pars, const_pars=const_pars, method=method)
}


bootstrap_tau <- function(X, const, timedep, indicators = c("Variance", "Autocorrelation", "Skew", "Kurtosis"), nboot=160, cpu=16, windowsize=round(length(X)/2), method=c("pearson", "kendall", "spearman")){
# Tau approach comparison
	taus <- lapply(indicators, function(stat){ 	tau_dist_montecarlo(X, const, timedep, signal=stat, nboot=nboot, cpu=cpu, method=method) })
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
	par(mfrow=c(m,n), oma=c(8,8,8,4), mar=c(0,0,0,0))

	for(j in 1:m){
		for(i in 1:n){
			if(j == m){ xaxt <- "s"
			} else { xaxt <- "n" }
			if(i > 1){ yaxt <- "n" 
			} else { yaxt <- "s" }

			plot(taus[[i]][[j]], show_p = show_p, show_error=show_error, xaxt = xaxt, yaxt=yaxt, ...)

			if(j==1) mtext(data_names[i], NORTH<-3, cex=par()$cex.lab, line=2) ## data labels on top row
			if(i==1){
				mtext(taus[[i]][[j]]$signal, WEST<-2, line=4) ## statistic name on first column
				mtext(expression(paste("Prob Density of ", tau)),
						WEST<-2, line=2, cex=.7*par()$cex.lab) ## statistic name on first column
			}	
			if(j==m & i==2) mtext(expression(paste(tau)), SOUTH<-1, line=4, cex=par()$cex.lab) ## x-axis label
		}
	}
}

