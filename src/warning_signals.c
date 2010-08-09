/**
 * @file warning_signals.c
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 * This implements one-dimensional (unstructured) stochastic population
 * dynamics for common functions, including those containing saddle node
 * bifurcations.  This file is still in early development, but will be used
 * to explore statistics of regime shifts.  
 * 
 * @section LICENCE
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */
#include "warning_signals.h"

/** Storage of statistics in window wraps around, so that once 
 * the sample window has been filled, it starts overwriting 
 * from the beginning.  This saves space but breaks the temporal
 * ordering of the data, which throws off some statistics (such as
 * autocorrelation or skew).  This function reorders the window so
 * that events are listed in temporal order.  Eventually might be 
 * best to store all data and avoid the extra computations!  */
void reorder(double * data, const int i, const size_t windowsize){
	int j,k;
	double tmp[windowsize];
	for(j=0;j<windowsize;j++){
		k = ( i + j ) % windowsize;
		tmp[j] = data[k];
	}
	for(j=0;j<windowsize;j++)
		data[j] = tmp[j];
}

/** Execute all tasks that occur on fixed interval schedule 
 * Currently this function isn't very api like, as it is 
 * closely tied to the logical structure of the record data structure 
 * which itself could be abstracted more */
void warning_fixed_interval(const double t, void * mypars, void * myrecord, int rep)
{
	pars * my_pars = (pars *) mypars;
	record * my_record = (record *) myrecord;
	size_t windowsize = my_record->windowsize;
	double * checkpts = my_pars->checkpts;

	/* Sample (print) system state at regular intervals */
	if(t > checkpts[0]){
		double * hist = my_record->hist;
		hist[my_pars->hist_index] = my_pars->n;


		if(t > my_record->sampletime){ /*  don't start printing until we've filled out the vector */
			int i = my_pars->time_index; //replicate 
			my_record->means[i] +=  gsl_stats_mean( hist, 1, windowsize);
			my_record->vars[i] +=  gsl_stats_variance( hist, 1, windowsize);
			reorder(hist, my_pars->hist_index, windowsize);
			my_record->skews[i] += gsl_stats_skew( hist, 1, windowsize);
			my_record->ar1[i] +=  gsl_stats_lag1_autocorrelation( hist, 1, windowsize);
			/* Lag N/2 autocorrelation found by splitting data in half */			
			my_record->arN[i] += gsl_stats_correlation(hist, 1, &(hist[windowsize/2]), 1, windowsize/2);

			my_record->a[i] = my_pars->a;
			++(my_pars->time_index);
		}
		my_pars->hist_index = (1+my_pars->hist_index) % windowsize;
		checkpts[0] += my_record->samplefreq;
	}

	/* Slowly change the bifurcation parameter */
	if(t > checkpts[1]){
		checkpts[1] += my_pars->pollute_rate;
		my_pars->a += my_pars->pollute_increment;	
	}


}

/** A function to reset pars to the initial conditions 
 * for ensemble runs.  Creates a copy so the original
 * initial conditions aren't destroyed while each sim
 * has its own copy to modify as system evolves.  */
void * warning_reset(const void * in)
{
	const pars * mypars = (const pars *) in;
	pars * my_pars = (pars *) malloc(sizeof(pars));
	my_pars->n = mypars->n;
	my_pars->K = mypars->K;
	my_pars->e = mypars->e;
	my_pars->a = mypars->a;
	my_pars->h = mypars->h;
	my_pars->time_index = mypars->time_index;
	my_pars->hist_index = mypars->hist_index;
	int i;
	for(i=0; i< N_EVENT_TYPES; i++){
		my_pars->checkpts[i] =	mypars->checkpts[i];
	}
	my_pars->start_polluting = mypars->start_polluting;
	my_pars->pollute_rate = mypars->pollute_rate;
	my_pars->pollute_increment = mypars->pollute_increment;
	return my_pars;
}

/** Stitching it all together */
void warning_signals(
	double * time, 
	double * a, 
	double * means, 
	double * vars,
	double * skews,
	double * ar1, 
	double * arN,	
	int * sample_time, ///< time window over which to compute stats
	double * sample_freq, ///< sampling frequency
	int * max_time, 
	int * n_ensembles, ///< number of ensembles
	double * start_polluting,
	double * pollute_rate,
	double * pollute_increment )
{
	/** Create a list of all event functions and their associated outcomes 
	 *  Order doesn't matter, but make sure outcomes are paired with the 
	 *  appropriate function! */
	event_fn rate_fn[N_EVENT_TYPES];	
	event_fn outcome[N_EVENT_TYPES];
	rate_fn[0] = &birth;
	outcome[0] = &birth_outcome;
	rate_fn[1] = &death;
	outcome[1] = &death_outcome;

	RESET reset_fn = &warning_reset;
	FIXED fixed_interval_fn = &warning_fixed_interval;


	record * my_record = ws_record_alloc(*sample_time, *sample_freq, *max_time);
	/** Allocate and intialize the parameters structure for the functions */
	int No = 572;
	int K = 1000;
	double e = 0.5;
	int Ao = 160;
	double h = 200;

	pars * my_pars = pars_alloc(No, K, e, Ao, h, 
								*start_polluting, 
								*pollute_rate, 
								*pollute_increment);
	size_t ensembles = *n_ensembles;

	gillespie(rate_fn, outcome, N_EVENT_TYPES, 
				my_pars, my_record, *max_time, 
				ensembles, reset_fn, fixed_interval_fn);

//	euler(my_pars, *max_time, stderr);
//	gslode(my_pars, *max_time, theory);
	
	double t = my_record->sampletime;
	int i=0;
//	fprintf(stderr, "\n");
	while(t<my_record->maxtime){
		t+=my_record->samplefreq;
		time[i] = t;
		a[i] = my_record->a[i];
		means[i] = my_record->means[i]/ensembles;
		vars[i] = my_record->vars[i]/ensembles;
		skews[i] = my_record->skews[i]/ensembles;
		ar1[i] = my_record->ar1[i]/ensembles;
		arN[i] = my_record->arN[i]/ensembles;
//		ens_mean[i] = my_record->ar1[i]/ensembles;
//		ens_var[i] = my_record->arN[i]/ensembles -gsl_pow_2( my_record->ar1[i]/ensembles );
/*		printf("%g %g %g %g %g %g %g\n",
				t,
				my_record->a[i],
				my_record->means[i]/ensembles, 
				my_record->vars[i]/ensembles,
				my_record->skews[i]/ensembles,
				my_record->ar1[i]/ensembles, //being used for ensemble mean
				my_record->arN[i]/ensembles -gsl_pow_2( my_record->ar1[i]/ensembles ) //being used for esemble variance
//				my_record->arN[i]/ensembles 
		);
*/
		++i;
	}

	free(my_pars);
	ws_record_free(my_record);
}

/* ws "main" function */
int ws(void)
{
	const size_t S = 500;
	double  time[S], a[S], means[S], vars[S], skews[S], ar1[S], arN[S];
	double sample_freq = 1;
	double start_polluting = 450;
	double pollute_rate = 1; 
	double pollute_increment = 1;
	int sample_time = 1; 
	int max_time = 500;
	int n_ensembles = 1;
	warning_signals(time, a, means, vars, skews, ar1, arN, &sample_time, 
					&sample_freq, &max_time, &n_ensembles, &start_polluting,
					&pollute_rate, &pollute_increment);
	int i;
	for(i=0; i<S; i++)
		printf("%g %g %g\n", time[i], a[i], means[i] );
	
	return 0;
}

