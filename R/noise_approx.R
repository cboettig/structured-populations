# A one-dimensional test case: b(n) = n, d(n) = n^2/K

require(odesolve)



oneD <- function(t, y, p){
	# \dot x = \alpha_1(x)
	yd1 <- p["r"]*y[1] *(1 - y[1]/p["K"])
	# \dot \sigma^2 = - \partial_x \alpha_1(x) \sigma^2 + \alpha_2(x)
	yd2 <- 2*p["r"]*(1-2*y[1]/p["K"] ) * y[2] + 
		p["r"]*y[1] *(1 + y[1]/p["K"])

	list(c(yd1, yd2))
}

logistic_parameters <- c(K=100, r=1)
times <- seq(0,10,by=.1)
yo <- c(10, 0)
o <- lsoda(yo, times, oneD, logistic_parameters)


# two-dimenstional test case: 
twoD <- function(t, y, p){


	# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"] * y[1] * (p["K"] - y[1] - y[2]) - p["d1"] * y[1] + p["c1"] * y[1] * y[2]
	# \dot y = \beta_1(x,y)
	yd2 <- p["b2"] * y[2] * (p["K"] - y[1] - y[2]) - p["d2"] * y[2] - p["c2"] * y[1] * y[2]
	# sigma_x^2
	yd3 <-  2*( p["b1"] * (p["K"] - 2*y[1] - y[2]) - p["d1"] + p["c1"] * y[2] ) * y[3]	+		# 2 \partial_x \alpha_1(x,y)  * \sigma_x^2
		(-p["b1"] * y[1]  + p["c1"] * y[1]  ) *  y[5]									+		# \partial_y \alpha_1(x,y)  * cov(x,y)
		p["b1"] * y[1] * (p["K"] - y[1] - y[2]) + p["d1"] * y[1] + p["c1"] * y[1] * y[2]		# \alpha_2
	# sigma_y^2
	yd4 <- 2*( p["b2"] * (p["K"] - y[1] - 2*y[2]) - p["d1"] + p["c2"] * y[1] ) * y[4]	+		# 2 \partial_y beta_1(x,y)  * \sigma_y^2
		(-p["b2"] * y[2]  + p["c2"] * y[2]  ) * y[5]									+		# \partial_x \beta_1(x,y) * cov(x,y)
		p["b1"] * y[1] * (p["K"] - y[1] - y[2]) + p["d1"] * y[1] + p["c1"] * y[1] * y[2]		# \beta_2
	# cov 
	yd5 <- (  2*( p["b1"] * (p["K"] - 2*y[1] - y[2]) - p["d1"] + p["c1"] * y[2] )  +			# (\partial_x \alpha_1 +
		( p["b2"] * (p["K"] - y[1] - 2*y[2]) - p["d1"] + p["c2"] * y[1] ) ) * y[5] +			#   \partial_y \beta_1 )* cov 
		(-p["b2"] * y[2]  + p["c2"] * y[2]  ) * y[3]	+										# \partial_x \beta_1 * sigma_x^2
		(-p["b1"] * y[1]  + p["c1"] * y[1]  ) * y[4]

	list(c(yd1, yd2, yd3, yd4, yd5))
}
k<-1000
#crowley_parameters <- c(b1=.11/k, b2=.6/k, d1=0.1, d2=.1, c1=0.1/k, c2=4/k, K=k)
crowley_parameters <- c(b1=.11/k, b2=.6/k, d1=0.1, d2=.1, c1=0.1/k, c2=3.9/k, K=k)
times <- seq(0,100,length=1000)
yo <- c(50, 450, 0, 0, 0) # (xo, yo, sigma_xo sigma_yo, cov)
o2 <- lsoda(yo, times, twoD, crowley_parameters)

png("diverge.png")
par(mfrow=c(2,1) )
plot(o2[,1], o2[,2], type='l', col="darkblue", ylim=c(0,1100), ylab = "expected abundance", xlab= "time")
lines(o2[,1], o2[,3], col="darkgreen")
legend("topleft", c("competitor", "colonist"), col=c("darkblue", "darkgreen"), pch='l')
plot(o2[,1], sqrt(o2[,5]), type='l', ylab = "stdev", xlab="time")
lines(o2[,1], sqrt(o2[,4]), col="darkblue")
lines(o2[,1], sqrt(o2[,5]), col="darkgreen")
dev.off()

# solve the system of equations
eqns <- function(t, x, p){
## add macroscopic equations for the state

## Covariance equations
	n <- sqrt( length(x) )
	y <- matrix(x, n,n)
	yd <- diag(H(x,p), n,n)
	for(i in 1:n){
		for(j in 1:n){
			if(i == j){
				yd[i,i] = yd[i,i] + 2*J(y,p)[i,j]*y[i,j]
			} else {
				yd[i,j] <- ( J(y,p)[i,i]+ J(y,p)[j,j])*y[i,j] + J(y,p)[j,i] * y[i,i] +  J(y,p)[i,j] * y[j,j]
			}
		}
	}
	list(as.numeric(yd))
}





p = c(	b=5, ue= 0, ul = .001, up = 0, ua = .003, ae = 1/3.8, al = 1/(20.2-3.8), ap = 1/(25.5-20.2), cle = 0.01, cap = 0.004, cae = 0.01)
	

J <- function(y,p){
 j <- c(-p["ue"]-p["ae"],	-p["cle"],	0,			p["b"] - p["cae"],
	  p["ae"],			-p["ul"],	0,			0,
	  0,				p["al"],	-p["up"],	0,
	  0,				0,			p["up"],	-p["ua"] )
 t(matrix(j,4,4))
}

H <- function(y,p){
	c(p["ue"]+p["ae"]+p["cle"]+p["b"]+p["cae"], p["ae"]+p["ul"], p["al"]+p["up"], p["up"]+p["ua"] )
}

