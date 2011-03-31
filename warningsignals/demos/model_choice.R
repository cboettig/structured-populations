## model_choice.R

## load libraries
require(warningsignals)
require(socialR)
gitcommit()

## load datafiles
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts

nboot <- 2000
cpu <- 16

## fit models
CaCO3_LSN <- fit_models(CaCO3, "LSN")
deut3_LSN <- fit_models(deut3, "LSN")
CaCO3_LTC <- fit_models(CaCO3, "LTC")
deut3_LTC <- fit_models(deut3, "LTC")

CaCO3_modelchoice <- montecarlotest(CaCO3_LTC$timedep, CaCO3_LSN$timedep, nboot=nboot, cpu=cpu)

deut3_modelchoice <- montecarlotest(deut3_LTC$timedep, deut3_LSN$timedep, nboot=nboot, cpu=cpu)

social_plot(plot(CaCO3_modelchoice), file="CaCO3.png", tag=paste(tags, "MC", "modelchoice CaCO3")) 

social_plot(plot(deut3_modelchoice), file="deut3.png", tag=paste(tags, "MC", "modelchoice CaCO3")) 

