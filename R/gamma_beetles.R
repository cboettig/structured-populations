## Simulate beetle dynamics


k <-10
eggs <- 1:k	
larva <- (k+1):(2*k)
pupa <- (2*k+1):(3*k)
adults <- 3*k+1

## Beetle model macroscopic eqns
b_gamma <- function(t,y,p){
	o <- numeric(adults);
	o[1] <- p["b"]*y[adults]
	o
}

## Beetle model macroscopic eqns
d_gamma <- function(t,y,p){
	o <- numeric(adults)

	o[eggs] <-  (p["ue"] + p["cle"]*sum(y[larva]) + p["cae"]*y[adults])*y[eggs]  
	o[larva] <- p["ul"]*y[larva]
	o[pupa] <- (p["up"]+p["cap"]*y[adults])*y[pupa] 
	o[adults] <- p["ua"]*y[adults]
	o
}
## Jacobian of f
J_gamma <- function(t,y,p){
	o <- matrix(0, adults, adults)
	d <- numeric(adults)

	#diagonal
	d[eggs] <- -p["ue"] - p["ae"] -p["cle"]*sum(y[larva]) - p["cae"]*y[adults]
	d[larva] <- -p["ul"] - p["al"]
	d[pupa] <- -p["up"] - p["ap"] -p["cap"]*y[adults]
	d[adults] <- -p["ua"] 
	diag(o) <- d

	# lower diagonal,
	o[matrix(c(eggs+1,eggs),k)] <-  p["ae"]
	o[matrix(c(larva+1,larva),k)] <- p["al"]
	o[matrix(c(pupa+1,pupa),k)] <- p["ap"]

	# cannibalism
	o[matrix(c(eggs, eggs+k),k)] <- -p["cle"]*y[eggs]
	o[matrix(c(eggs, rep(adults,k)),k)] <- -p["cae"]*y[eggs]
	o[1, adults] <- o[1,adults] + p["b"]*y[adults]
	o[matrix(c(pupa, rep(adults,k)),k)] <- -p["cap"]*y[pupa]
	o
}




## Beetle model macroscopic eqns
T_gamma <- function(t,y,p){
	o <- matrix(0,adults,adults)
	o[matrix(c(eggs,eggs+1),  k)] <-  p["ae"]*y[eggs]
	o[matrix(c(larva,larva+1),k)] <- p["al"]*y[larva]
	o[matrix(c(pupa,pupa+1), k)] <- p["ap"]*y[pupa]
	o
}

collapse_gamma_classes <- function(beetle_data, k){
	# Collapse the pseudo-classes
	data <- matrix(0, length(times), 9)

	#time
	data[,1] = beetle_data[,1]

	# means
	for( i in 1:3) {
		data[,i+1] = rowSums(beetle_data[,(2+k*(i-1)):(1+k*i)] )
	}
	data[,5] = beetle_data[,3*k+2]
	vars_start <- 3*k+3
	for( i in 1:3) {
		data[,5+i] = rowSums(beetle_data[,(vars_start+k*(i-1)):(vars_start+1+k*i)] )
	}
	data[,9] = beetle_data[,1+2*(3*k+1)]
	data
}




