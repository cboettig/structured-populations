## FIGURE 2
source("load_CaCO3.R")
CaCO3 <- X
source("load_deut.R")
deut3 <- data[[3]]$X_ts

CaCO3_m <- fit_models(CaCO3, "LSN")
CaCO3_taus <- bootstrap_tau(CaCO3_m$X, CaCO3_m$const, CaCO3_m$timedep)
deut3_m <- fit_models(deut3, "LSN")
deut3_taus <- bootstrap_tau(deut3_m$X, deut3_m$const, deut3_m$timedep)

save(list=ls(), file="real_data_examples.Rdat")

social_plot(
	plot.bootstrap_tau(data.frame(CaCO3_taus, deut3_taus)),
	tags = "warningsingals stochpop"
	)


