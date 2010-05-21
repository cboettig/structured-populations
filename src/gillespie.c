/**
 * @file gillespie.c
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 * To use this function, one must write functions for determining 
 * the rates of all stochastic events, and for determining the outcomes of 
 * each of those events.  Usually the user will also specify events that 
 * happen on a fixed (deterministic) interval schedule, such as sampling 
 * the state of the system.  This is done by specifying the fixed_interval_tasks
 * function.  
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


#include "gillespie.h"



/** @short Simulates a continuous time stochastic process exactly.
 *
 * The gillespie function computes the rates of all possible events
 * for the current system state.  The next event occurs after a time
 * interval t which is an exponential random variable with rate equal
 * to the sum of all these rates.  The event  which occurs at that
 * time is drawn at random from the events list weighted by the rates.  
 * More presicely, if RATES is an array with the cumulative sum of all
 * weights divided by the sum of all rates, then the chosen event is
 * the one corresponding to the first number in the cumulative sum 
 * larger than a uniform random variable (on 0,1). 
 *
 * @param rate_fn 
 *			An array of function pointers for calculating the 
 *			rates at which each possible event can occur.  Order is irrelevant but 
 *			must match the order used in outcomes array. Rate functions must 
 *			take a single argument, void pointer promoted to necessary type
 * @param outcome 
 *			An array of fn ptrs for creating the outcome of each event. 
 *			Order must match rate_fn array and take the same argument.
 * @param n_event_types 
 *			Number of different types of events
 * @param my_pars 
 *			A pointer to a structure of type pars for specifying 
 *			the parameters and variables. This is the argument to rates and 
 *			outcomes fns.  my_pars should also contain an array of times at 
 *			which fixed events such as sampling events occur.   
 * @param max_time 
 *			Length of time to simulate each replicate
 * @param ensembles 
 *			Number of replicates to simulate
 *
 * @return void.  Usually system status is printed out by the fixed_interval_tasks fn.
 * */
void 
gillespie(	const event_fn * rate_fn,
			const event_fn * outcome,
			const size_t n_event_types,
			void * inits,	
			void * my_record,
			const size_t max_time,	
			const size_t ensembles,
			RESET reset_fn,
			FIXED fixed_interval_fn
		)
{

	/* intialize random number generator and seed with timestamp */
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);
	#ifdef _RNDSEED
	gsl_rng_set(rng, time(NULL));
	#endif

	/* Gillespie simulation variables */
	int l,i,check;
	double lambda, t, tmp;
	double * rates_data;

    /* private variables are accessible inside the parallel region only
	 * Dynamically allocated private arrays must be declared inside 
	 * the parallel region.  Currently this uses the pars_alloc function 
	 * to allocate space in user-defined way without assuming the function type */
	#pragma omp parallel shared(rng, rate_fn, outcome, my_record, inits, reset_fn, fixed_interval_fn) \
	private(lambda, t, tmp, i, check, rates_data,l)
	{
		/* The vector to store cumulative sum of rates */
		rates_data = (double *) calloc (n_event_types,sizeof(double) );

		/* Loop over ensembles, will be parallelized if compiled with -fopenmp */
		#pragma omp for
		for(l = 0; l < ensembles; l++){
			void * my_pars = reset_fn(inits);
			t = 0;
			check=0;
			while(t < max_time){
				/* calculate time until next event */
				rates_data[0] = rate_fn[0](my_pars);
				for(i=1;i<n_event_types;i++){
					rates_data[i] = rates_data[i-1] + rate_fn[i](my_pars);
				}
				lambda = rates_data[n_event_types-1]; 
				t += gsl_ran_exponential(rng, 1/lambda);
			
				/* Events such as sampling occur at regular intervals */
				fixed_interval_fn(t, my_pars, my_record, l);

				/* Determine if event is a birth or death */
				tmp = gsl_rng_uniform(rng);
				for(i=0;i<n_event_types;i++){
					if( tmp < rates_data[i]/lambda ){
						check = outcome[i](my_pars);
						break;
					} 
				}
			if(check) break;

		   }/* end evolution */
			free(my_pars);
		}/* end ensembles */
		free(rates_data);
	} /* end parallel */
	gsl_rng_free(rng);
}



