source("load_CaCO3.R")
social_plot(all_indicators(X), tags=tags, file="CaCO3.png")

source("load_deut.R")
## Glaciation I, processed
social_plot(all_indicators(data[[1]]$X_ts), tags=tags, file="DeutI.png")
## Glaciation II, processed
social_plot(all_indicators(data[[2]]$X_ts), tags=tags, file="DeutII.png")
## Glaciation III, processed
social_plot(all_indicators(data[[3]]$X_ts), tags=tags, file="DeutIII.png")
## Glaciation IV, processed
social_plot(all_indicators(data[[4]]$X_ts), tags=tags, file="DeutIV.png")


source("load_deut.R")
## Glaciation I, unprocessed
social_plot(all_indicators(data[[1]]$raw_ts), tags=tags)
## Glaciation II, unprocessed
social_plot(all_indicators(data[[2]]$raw_ts), tags=tags)
## Glaciation III, unprocessed
social_plot(all_indicators(data[[3]]$raw_ts), tags=tags)
## Glaciation IV, unprocessed
social_plot(all_indicators(data[[4]]$raw_ts), tags=tags)






##### should be able to pass raw time series as such, but no luck...
## Glaciation I, unprocessed
social_plot(all_indicators(data[[1]]$X), tags=tags)
## Glaciation II, unprocessed
social_plot(all_indicators(data[[2]]$X), tags=tags)
## Glaciation III, unprocessed
social_plot(all_indicators(data[[3]]$X), tags=tags)
## Glaciation IV, unprocessed
social_plot(all_indicators(data[[4]]$X), tags=tags)


