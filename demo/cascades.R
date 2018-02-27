# cascades.R
# dN/dt = f(N,P)
# dP/dt = g(N,P)

require(deSolve)

parms <- c(r=1, K=500, alpha=0.01, beta=0.25, mu=0.2, sigma=.001)


env <- function(t, parms){
  if(t>50 & t<100) 
    parms["K"] <- parms["K"] - 4*t
  if(t>100) 
    parms["K"] <- parms["K"] - 400
  parms
}

system <- function(env, i){
  f <- function(t, y, parms){
    parms <- env(t, parms)
    y1_p <- parms['r']*y[1]*(1-y[1]/parms["K"]) - parms["alpha"]*y[1]*y[2] #+ parms["sigma"]*rnorm(1)*y[1]
    y2_p <- parms["alpha"]*parms["beta"]*y[1]*y[2] - parms["mu"]*y[2]
    # Note: second return element is for global parms as funct of time t
    list(c(y1_p, y2_p), NULL) 
  }

  Nhat <- parms["mu"]/(parms["alpha"]*parms["beta"])
  Phat <- parms["r"]/parms["alpha"] * (1- Nhat/parms["K"])
  print(paste("Initial N_hat =", Nhat, "P_hat =", Phat))
  y <- c(N=Nhat, P=Phat)
  times <- seq(0,200,length.out=100)
  out <- lsoda(y, times, f, parms=parms,  rtol=1e-4, atol=1e-4)
  enviro <- sapply(times, env, parms)

  par(mfrow=c(2,1))
  plot(out[,1], out[,2], type="l", col="green", 
       ylim = c(min(out[,2], out[,3]), max(out[,2], out[,3])), lwd=2)
  lines(out[,1], out[,3], col="red", lwd=2)
  legend("topright", c("prey", "predator"), col=c("green", "red"), lty=1)
  plot(times, enviro[i,], col="black", lwd=2, type="l")

  parms <- enviro[,100]
  Nhat <- parms["mu"]/(parms["alpha"]*parms["beta"])
  Phat <- parms["r"]/parms["alpha"] * (1- Nhat/parms["K"])
  print(paste("Final N_hat =", Nhat, "P_hat =", Phat))
}

png("K.png")
system(env, 2)
dev.off()

env <- function(t, parms){
  if(t>50 & t<100) 
    parms["mu"] <- parms["mu"] + t/2000
  if(t>100) 
    parms["mu"] <- parms["mu"] + 100/2000
  parms
}

png("mu.png")
system(env,5)
dev.off()

require(socialR)
gitopts <- c(user = "cboettig", repository = "structured-populations", dir = "demo", raw = FALSE, diff = FALSE)
upload("mu.png", script="cascades.R", gitaddr=gitcommit("cascades.R", gitopts, msg="autocommit on figure upload"))
upload("K.png", script="cascades.R", gitaddr=gitcommit("cascades.R", gitopts, msg="autocommit on figure upload"))

