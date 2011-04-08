# Simulate a dataset from the full individual, nonlinear model
T<- 1000
n_pts <- 40
require(warningsignals)
pars = c(Xo = 730, e = 0.5, a = 100, K = 1000, h = 200, 
    i = 0, Da = .09, Dt = 0, p = 2)
sn <- saddle_node_ibm(pars, times=seq(0,T, length=n_pts))
ibm_critical <- ts(sn$x1,start=sn$time[1], deltat=sn$time[2]-sn$time[1])

plot(ibm_critical, cex.lab=2, cex.axis=2, lwd=1, xlab="time", ylab="pop")

# Stable
pars = c(Xo = 730, e = 0.5, a = 100, K = 1000, h = 200, 
    i = 0, Da = 0, Dt = 0, p = 2)
sn <- saddle_node_ibm(pars, times=seq(0,T, length=n_pts))
ibm_stable  <- ts(sn$x1,start=sn$time[1], deltat=sn$time[2]-sn$time[1])

plt <- function(){
plot(ibm_critical, cex.lab=2, cex.axis=2, lwd=1, xlab="time", ylab="pop")
lines(ibm_stable)
}

require(socialR)
social_plot(plt(), tags="ibm stochpop warningsignals")

save(list=c("ibm_critical", "ibm_stable"), file="ibm_sims.Rdat")

