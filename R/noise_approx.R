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
		(-p["b1"] * y[1]  + p["c1"] * y[1]  ) * y[4]											# \partial_y \alpha_1 * sigma_y^2

	list(c(yd1, yd2, yd3, yd4, yd5))
}






linnoise <- function(t,y,p, f, g, J){
	n <- length(y)				# d states, d variances, (d^2-d)/2 covariances
	d <- -3/2 + sqrt(9/4+2*n)	# dimension of the system
	yd <- numeric(n)

	jac <- J(t,y,p) 

	# macroscopics
	yd[1:d] <- f(t,y,p)

	# variances
	yd[ (d+1):(2*d) ] <- 2 * jac %*% y[ (d+1):(2*d) ] + g(t,y,p)
	
	# covariances
	k <- 2*d+1
	for (i in 1:(d-1)){
		for (j in (i+1):d){
				# Computes ( J(i,i) + J(j,j) )*Cov(i,j)  + J(j,i)*Cov(i,i) + J(i,j) Cov(j,j)
				yd[k] <- ( jac[i,i]+ jac[j,j])*y[k] + jac[j,i] * y[i+d] +  jac[i,j] * y[j+d]
				k = k+1
		}
	}
	list(yd);
}


f_crowley <- function(t,y,p){
	# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"] * y[1] * (p["K"] - y[1] - y[2]) - p["d1"] * y[1] + p["c1"] * y[1] * y[2]
	# \dot y = \beta_1(x,y)
	yd2 <- p["b2"] * y[2] * (p["K"] - y[1] - y[2]) - p["d2"] * y[2] - p["c2"] * y[1] * y[2]
	c(yd1, yd2)
}

g_crowley <- function(t,y,p){
# \dot x = \alpha_1(x,y)
	yd1 <- p["b1"] * y[1] * (p["K"] - y[1] - y[2]) + p["d1"] * y[1] + p["c1"] * y[1] * y[2]
	# \dot y = \beta_1(x,y)
	yd2 <- p["b2"] * y[2] * (p["K"] - y[1] - y[2]) + p["d2"] * y[2] + p["c2"] * y[1] * y[2]
	c(yd1, yd2)
}


J_crowley <- function(t,y,p){
	j <- c(		( p["b1"] * (p["K"] - 2*y[1] - y[2]) - p["d1"] + p["c1"] * y[2] ),	(-p["b1"] * y[1]  + p["c1"] * y[1]  ), 
			(-p["b2"] * y[2]  + p["c2"] * y[2]  ),								( p["b2"] * (p["K"] - y[1] - 2*y[2]) - p["d1"] + p["c2"] * y[1] ) ) 
	t(matrix(j,2,2))
}

# lsoda needs a function of only t,y,p


k<-1000
crowley_parameters <- c(b1=.11/k, b2=.6/k, d1=0.1, d2=.1, c1=0.1/k, c2=3.9/k, K=k)
times <- seq(0,10,length=100)
yo <- c(50, 450, 0, 0, 0) # (xo, yo, sigma_xo sigma_yo, cov)


o2 <- lsoda(yo, times, twoD, crowley_parameters)


eqns <- function(t,y,p){ linnoise(t,y,p, f_crowley, g_crowley, J_crowley) }
testcase <- lsoda(yo, times, eqns, crowley_parameters)









# y is given as  c( state, variances, covariances )  states = d, then length = d+d+(d^2-d)/2 = 1.5*d+d^2/2 = n
#d <- 4
#n <- 1.5*d+d^2/2

#Yo <- numeric(n)
#Yo[1] <- 100


## Beetle model parameters
p = c(	b=5, ue= 0, ul = .001, up = 0, ua = .003, ae = 1/3.8, al = 1/(20.2-3.8), ap = 1/(25.5-20.2), cle = 0.01, cap = 0.004, cae = 0.01)

## Beetle model macroscopic eqns
f <- function(t,y,p){
	f1 <- p["b"]*y[4]  - (p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4] )*y[1]  
	f2 <- p["ae"]*y[1] - (p["ul"] + p["al"])*y[2]
	f3 <- p["al"]*y[2] - (p["up"] + p["ap"])*y[3] 
	f4 <- p["ap"]*y[3] - p["ua"]*y[4] 
	c(f1, f2, f3, f4)
}

## Beetles second jump moment
g <- function(t,y,p){
	g1 <- p["b"]*y[4]  + (p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4] )*y[1]  
	g2 <- p["ae"]*y[1] + (p["ul"] + p["al"])*y[2]
	g3 <- p["al"]*y[2] + (p["up"] + p["ap"])*y[3] 
	g4 <- p["ap"]*y[3] + p["ua"]*y[4] 
	c(g1, g2, g3, g4)
}

## Jacobian of f
J <- function(t,y,p){
 j <- c( -(p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4]),	-p["cle"]*y[1],	0,			p["b"] - p["cae"]*y[1],
	  p["ae"],			-p["ul"]-p["al"],	0,			0,
	  0,				p["al"],	-p["up"]-p["ap"],	0,
	  0,				0,			p["ap"],	-p["ua"] )
 t(matrix(j,4,4))
}


	

