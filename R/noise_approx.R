# solve the system of equations
eqns <- function(t, x, p){
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

