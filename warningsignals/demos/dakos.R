# data from: 
deut <- read.table("deutnat.txt")

# divisions from: Dakos et. al.
g1 <- 1:591
g2 <- 592:(259+592)
g3 <- (1+259+592):(1+259+592+150)
g4 <- (2+259+592+150):(2+259+592+150+127)


glaciationI <- c(deut$V2[g1], deut$V3[g1])
glaciationII <- c(deut$V2[g2], deut$V3[g2])
glaciationIII <- c(deut$V2[g3], deut$V3[g3])
glaciationIV <- c(deut$V2[g4], deut$V3[g4])


caco3 <- read.table("caco3.txt")
caco3 <- data.frame("MYYrs"=-caco3$V1, "CaCO3"=caco3$V2)
dat <- caco3$MYYrs > -40 & caco3$MYYrs < -34
myyrs <- caco3$MYYrs[dat]
end=myyrs[1]; start=myyrs[length(myyrs)]


## Incorrect time intervals
X <- ts(caco3[dat,2],start=start,end=end, freq=length(dat)/(end-start))


plot(caco3[dat,1], caco3[dat,2])



