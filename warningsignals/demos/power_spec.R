#power_spec.R
require(sde)





n<- 2^10
dt<-1/2^2
t <- seq(0,dt*(n-1), by=dt)
X <- sin(5*t) + rnorm(n, .01)
pwr <- abs(fft(dt*X))^2/(2*pi*dt*n)
nyquist<-1/(2*dt)  # half the sampling frequency
w <- nyquist*seq(0, 2*pi, length=n/2)
plot(w, pwr[1:(n/2)])
abline(v=5)




a1p <- 1
a2 <- 1

n <- 2^10
d <- expression( -1*x )
s <- expression( 1 )
dt <- 1
T <- n*dt
reps<- 100

X <- sde.sim(drift = d, sigma = s, N=(n-1), delta=dt, M=reps)
X <- matrix(X, nrow=n)
pwr<- sapply(1:reps, function(i){
	y <- abs( fft(dt*X[,i]) )^2/ (2*pi*dt*n)
	y[1:(n/2)]
})
nyquist=.5/dt
freq <- nyquist*seq(0,2*pi, length=n/2)
plot(freq, rowMeans(pwr))
lorentz <- function(w, alpha2, alpha1_p){ alpha2/((alpha1_p^2+w^2)) }
lines(freq, lorentz(freq, a2, a1p), col="green")


lines(freq, lorentzian(freq, a2, a1p), col="blue")




# explicit conjugate transpose averaging, with or without replicates
X_tilde <- matrix(fft(X),nrow=reps)*dt  ## integrate times dt
transform <- sapply(1:n, function(i) t(Conj(X_tilde[,i])) %*% X_tilde[,i]  ) / reps / (dt*n)  
power <- c(transform[(n/2+1):n], transform[1:(n/2)])
w_p <- seq(-pi/dt, pi/dt, length=n)
plot(w_p, power)
lines(w_p, pi*lorentzian(w_p, a2, a1p), col="blue", lwd=3, type="l", lty=2 )




# no-replicate norm square
transform <- abs(fft(dt*X)^2) / (2*pi*dt*n)
power <- c(transform[(n/2+1):n], transform[1:(n/2)])
w_p <- seq(-pi/dt, pi/dt, length=n)
plot(w_p, power)

### Shows that the FFT of the correlation function matches the ##
### Lorentzian function in the frequency domain				   ##

# Define the Correlation function 
C <- function(tau, alpha2, alpha1_p) alpha2/abs(2*alpha1_p)*exp(-abs(alpha1_p)*abs(tau))
# Define the time domain
time <- seq(-dt*n, dt*n, length=2*n)

# Fourrier Transform the Correlation function and reshape
S_from_C <- abs(fft(C(time, a2, a1p)))
end <- length(S_from_C)
S_from_C <- c(S_from_C[(1+end/2):end], S_from_C[0:(end/2)] )

# define the Lorentzian 
lorentzian <- function(w, alpha2, alpha1_p){ alpha1_p/(pi*(alpha1_p^2+w^2)) }
# Define the frequency domain
w = seq(-pi/dt, pi/dt, length=2*n)

# Plot the Lorentzian and fft(C) 

plot(w_p, lorentzian(w_p, a2, a1p), col="blue", lwd=3, type="l", lty=2 )
lines(w,dt*S_from_C)
points(w_p, power)

l<-lorentzian(w_p, a2, a1p)
test<-X_tilde/sqrt(2*pi*dt*n*sqrt(l))
qqplot(test)

##### Try example with simpler simulations: AR(1) in (x), Euler version of OU in (y)

n <- 2^12
gamma <- 8
D <- .3
dt <- 0.05
S <- function(w) 2*D/((gamma)^2+w^2)
## autocorr example
x <- numeric(n)
y<-x

pwr <- sapply(1:10, function(j) {
for(i in 2:n){
#	x[i] <- x[i-1]*.5 +sqrt(1-.5^2)* rnorm(1)
#	y[i] <- y[i-1]*(1-gamma*dt) +sqrt(2*D*dt)*rnorm(1)
	y[i] <- y[i-1]*exp(-gamma*dt) + rnorm(1, mean=0, sd=sqrt( (D/gamma)*(1-exp(-2*gamma*dt)))) 
}
fy <- abs(fft(dt*y))^2/(2*pi*n*dt)
fy <- c(fy[(n/2+1):n], fy[1:(n/2)])
fy
})

fy <- rowMeans(pwr)
#  is nyquist freq, since sampled every unit t.  also means n = T
nyquist <- 1/(2*dt)
omega <- seq(-2*pi*nyquist, 2*pi*nyquist, length=length(fy))
plot(omega, fy)
lines(omega, S(omega) )





fx <- abs(fft(x))^2/(2*pi*n*dt)
fx <- c(fx[(n/2+1):n], fx[1:(n/2)])
#  is nyquist freq, since sampled every unit t.  also means n = T
nyquist <- dt/2
omega <- seq(-2*pi*nyquist, 2*pi*nyquist, length=length(fy))
plot(omega, fx)

lines(lorentzian(omega, 1, .5))




sigma <- a2/(2*a1p)
fit_pars <- c(a2=a2, a1p=a1p, sigma=sigma)

f2 <- lorentzian(w_p, fit_pars["a2"], fit_pars["a1p"] )
plot(power-f2)
plot(power-f2 - fit_pars["sigma"]^2)

f<- sqrt(f2)

plot( (power-f2-fit_pars["sigma"]^2)/(2*fit_pars["sigma"]*f) )

plot(f+rnorm(length(f), sd=sigma) )


sigma <- 1
fit_pars <- c(a2=a2, a1p=a1p, sigma=sigma)
power_pdf <- function(w, fit_pars){
	f <- sqrt(lorentzian(w, fit_pars["a2"], fit_pars["a1p"]))
	f + dnorm(w, mean=0, sd=fit_pars["sigma"])
}

plot(w, power_pdf(w, fit_pars))










### Shows that the inverse FFT of the Lorentzian is the ##
### correlation function--- HMM... missing something	##


l <- lorentzian(w, a2, a1p)
end <- length(l)
C_from_S <- abs(fft(l, inverse=TRUE))
C_from_S <- c(C_from_S[(1+end/2):end], C_from_S[0:(end/2)] )
plot(time, C_from_S)
lines(time, C(time, a2, a1p)/dt^2)




#w = seq(-pi*end/(delta*n),pi*end/(2*delta*n), length=end);

#smooth_pwr <- loess.smooth(w_p, power, evaluation=length(w_p), family="gaussian", span=.05)
#lines(smooth_pwr)


corfn <- fft(power, inverse = T)






### with averaging 
m_samples <- 2^4
n_length <- n/m_samples  # number of points in the dataset
Y <- matrix(X@.Data, n_length, m_samples)
transform <- fft(Y)
tmp_power <- rowMeans( abs(transform)^2 )/n_length	# is this divided by n_length or m_samples??  is ths rowMeans or colMeans??
power <- c(tmp_power[(n_length/2+1):n_length], tmp_power[1:(n_length/2)])
w_p <- seq(-pi, pi, length=n_length)
plot(w_p, power)






