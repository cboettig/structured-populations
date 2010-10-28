#power_spec.R
require(sde)

a1p <- 150
a2 <- 2

n <- (2^10)
d <- expression( -150*x )
s <- expression( 2 )
dt <- .001
T <- n*dt

X <- sde.sim(drift = d, sigma = s, N=(n-1), delta=dt)
transform <- abs((fft(X@.Data))^2) /n^2
power <- c(transform[(n/2+1):n], transform[1:(n/2)])
w_p <- seq(-n*pi/T, n*pi/T, length=n)
plot(w_p, power)



### Shows that the FFT of the correlation function matches the ##
### Lorentzian function in the frequency domain				   ##

# Define the Correlation function 
C <- function(tau, alpha2, alpha1_p) alpha2/abs(2*alpha1_p)*exp(-abs(alpha1_p)*abs(tau))
# Define the time domain
time <- seq(-dt*n_length, dt*n_length, length=2*n_length)

# Fourrier Transform the Correlation function and reshape
S_from_C <- abs(fft(C(time, a2, a1p)))
end <- length(S_from_C)
S_from_C <- c(S_from_C[(1+end/2):end], S_from_C[0:(end/2)] )

# define the Lorentzian 
lorentzian <- function(w, alpha2, alpha1_p){ alpha2/(alpha1_p^2+w^2) }
# Define the frequency domain
w = seq(-pi/dt, pi/dt, length=2*n_length)

# Plot the Lorentzian and fft(C) 
plot(w, lorentzian(w, a2, a1p), col="blue", lwd=3, type="l", lty=2 )
lines(w,dt*S_from_C)

points(w_p, power)


### Shows that the inverse FFT of the Lorentzian is the ##
### correlation function--- HMM... missing something	##


l <- lorentzian(w, a2, a1p)
end <- length(l)
C_from_S <- abs(fft(l, inverse=TRUE))
C_from_S <- c(C_from_S[(1+end/2):end], C_from_S[0:(end/2)] )
plot(time, C_from_S)
lines(time, C(time, a2, a1p)/dt^2)




#w = seq(-pi*end/(delta*n_length),pi*end/(2*delta*n_length), length=end);

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






