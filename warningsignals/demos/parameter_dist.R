load("35563325713.Rdat")
require(socialR)
require(warningsignals)


plot_m <- function(pow, X=NULL){
    test <- pow$test_par_dist
    shade_plot <- function(par, lower=-Inf, upper=Inf){
# disgard nonsense values
        par <- par[par>lower & par < upper]
        dist <- density(par)
        plot(dist, type="n")
        polygon(dist$x, dist$y, col="gray", border="darkgray")
    }
    upper <- 1
    lower <- -1
    if(!is.null(X)){
        lower <- -10*length(X)/max(time(X))
        upper <- -lower
    }
    shade_plot(test[2,], lower, upper)
}


par(mfrow=c(1,2))
plot(density(deterior_mc$test_par_dist[2,], bw=1), xlim=c(-20,1))
plot(density(constant_mc$test_par_dist[2,c(-173,-1717)]))


## parameter distributions
par_names <- c("Ro", "m", "theta", "sigma")
plot_parameters <- function(pow){
    test <- pow$test_par_dist
    null <- pow$null_par_dist

    shade_plot <- function(par, lower=-Inf, upper=Inf, ...){
        par <- par[par>lower & par < upper]
        dist <- density(par)
        plot(dist, type="n", ...)
        polygon(dist$x, dist$y, col="gray", border="darkgray")
    }

    n <- dim(test)[1]
    par(mfrow=c(2,n))
    for(i in 1:n ){
        shade_plot(test[i,], main=par_names[i])
    }
    m <- dim(null)[1]
    for(i in 1:m ){
        if(i==2) plot(test[2,], type="n")
        shade_plot(null[i,], main="")
    }
}

social_plot(plot_parameters(deterior_mc), tag="sim deterior warningsignals stochpop parameters")
social_plot(plot_parameters(constant_mc), tag="sim constant warningsignals stochpop parameters")
social_plot(plot_parameters(deut3_mc), tag="sim deut3 warningsignals stochpop parameters") 





load("35575554231.Rdat")
converge_plot <- function(pow, ...){ 
    converged <- (pow$test_par_dist[5,] == 0)
    dens <- density(pow$test_par_dist[2,converged])
    plot(dens, type="n", ...)
    polygon(dens$x, dens$y, col="gray", border="darkgray")
}

plt <- function(){
par(mfrow=c(1,2))
    converge_plot(deterior_mc, main = "Deteriorating",
                  xlab="estimated stablization loss rate, m")
    converge_plot(constant_mc, main = "Constant", 
                  xlab="estimated stablization loss rate, m")
}
social_plot(plt(), tags = "warningsingals stochpop", 
			height=480, width=480*2)


cairo_pdf(file="parameter_dist.pdf", width=7, height=7/2)
    plt()
dev.off()
