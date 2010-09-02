
## Beetle model macroscopic eqns
b_beetles_one <- function(t,y,p){
	f1 <- p["b"]*y[4] 
	f2 <- p["ae"]*y[1]
	f3 <- p["al"]*y[2]
	f4 <- p["ap"]*y[3]
	c(f1, f2, f3, f4)
}

## Beetle model macroscopic eqns
d_beetles_one <- function(t,y,p){
	f1 <-  (p["ae"] + p["ue"] + p["cle"]*y[2] + p["cae"]*y[4] )*y[1] 
	f2 <- (p["al"]+p["ul"])*y[2]
	f3 <- (p["ap"]+p["up"]+p["cap"]*y[4])*y[3] 
	f4 <- p["ua"]*y[4] 
	c(f1, f2, f3, f4)
}

## Jacobian of f
J_beetles_one <- function(t,y,p){
 j <- c(	-(p["ue"] + p["ae"] + p["cle"]*y[2] + p["cae"]*y[4]),
			-p["cle"]*y[1],
			0, 
			p["b"] - p["cae"]*y[1],

			p["ae"],
			-p["ul"]-p["al"], 
			0, 
			0,

			0, 
			p["al"], 
			-p["up"]-p["cap"]*y[4]-p["ap"], 
			-p["cap"]*y[3],

			0,
			0,
			p["ap"],
			-p["ua"] )
 t(matrix(j,4,4))
}




## Beetle model macroscopic eqns
T_beetles_one <- function(t,y,p){
	t(matrix(c(
		0,		0,				0,				0,
		0,		0,				0,				0,
		0,		0,				0,				0,
		0,		0,				0,				0
		), 4,4))
}



