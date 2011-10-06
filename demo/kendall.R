# example regarding statistics in Dakos et al 2008
require(stochPop)
with_crash <- warning_signals(sample_time = 50, max_time=300, n_ensembles=1, start_polluting=200, pollute_timestep = .1, pollute_increment=.05)
plot(with_crash)

no_crash <- warning_signals(sample_time = 50, max_time=300, n_ensembles=1, start_polluting=400, pollute_timestep = .1, pollute_increment=.05)
plot(no_crash)


