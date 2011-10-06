# cascades.R
# dN/dt = f(N,P)
# dP/dt = g(N,P)
rm(list=ls())
require(odesolve)

parms <- c(r=5, K=1500, alpha=0.2, beta=0.75, mu=0.02, H=1500, a=1)


env <- function(t, parms){
  if(t>500 & t<1000) 
    parms["K"] <- parms["K"] - t
  if(t>1000) 
    parms["K"] <- parms["K"] - 1000
  parms
}

# functional response
g <- function(y,parms){
  #y[1] #linear
  y[1]^parms["a"]/(parms["H"]^parms["a"]+y[1]^parms["a"]) #saturating
}

  f <- function(t, y, parms){
    parms <- env(t, parms)
    y1_p <- parms['r']*y[1]*(1-y[1]/parms["K"]) - parms["alpha"]*y[1]*y[2] 
    y2_p <- parms["alpha"]*parms["beta"]*y[2]*g(y,parms) - parms["mu"]*y[2]
    # Note: second return element is for global parms as funct of time t
    list(c(y1_p, y2_p), NULL) 
  }

  # \mu H + \mu y_1 = \alpha \beta y_1
  # \mu H = (\alpha \beta - mu) y_1
  # r (1- y_1/K) = \alpha y_2 


  Nhat <- function(parms) parms["mu"]*parms["H"]/ ( parms["alpha"]*parms["beta"]-parms["mu"] )
  Phat <- function(parms) parms["r"]/parms["alpha"] * (1- Nhat(parms)/parms["K"])
  print(paste("Initial N_hat =", Nhat(parms), "P_hat =", Phat(parms)))

  y <- c(N=Nhat(parms), P=Phat(parms))


  times <- seq(0,2000,length.out=100)
  out <- lsoda(y, times, f, parms=parms,  rtol=1e-4, atol=1e-4)
  enviro <- sapply(times, env, parms)

  par(mfrow=c(2,1))
  plot(out[,1], out[,2], type="l", col="green", 
       ylim = c(min(out[,2], out[,3]), max(out[,2], out[,3])), lwd=2)
  lines(out[,1], out[,3], col="red", lwd=2)
  legend("topright", c("prey", "predator"), col=c("green", "red"), lty=1)
  plot(times, enviro[2,], col="black", lwd=2, type="l")

  parms <- enviro[,100]
  print(paste("Final N_hat =", Nhat(parms), "P_hat =", Phat(parms)))


#png("K2.png")
#system(env, 2)
#dev.off()

#png("mu2.png")
#system(env,5)
#dev.off()
