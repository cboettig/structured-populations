load("5572028881.Rdat")
require(socialR)
require(warningsignals)
indicators <- c("Variance", "Autocorrelation", "Skew", "Kurtosis")
nboot <- 2
cpu <-2 
ibm_tags="warningsignals stochpop ibm" 

deterior_mc <- 
		montecarlotest(deterior_m$const, deterior_m$timedep, 
		cpu=cpu, nboot=nboot, GetParNames=FALSE)


