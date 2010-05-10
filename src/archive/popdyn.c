/**
 * @file popdyn.c
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


/** Specify all event functions.  Functions should return the rate as type double
 * and take the sole argument as a pointer to the pars structure.  */
double death(void * mypars)
{ 
	pars * my_pars = (pars *) mypars;
	return my_pars->n * my_pars->e + my_pars->a; 
}
double birth(void * mypars)
{ 
	pars * my_pars = (pars *) mypars;
	return my_pars->e * my_pars->K * gsl_pow_2( my_pars->n ) / 
			(gsl_pow_2(my_pars->n) + gsl_pow_2(my_pars->h) ); 
}



/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double death_outcome(void * mypars)
{ 
	pars * my_pars = (pars *) mypars;
	--(my_pars->n);
	if(my_pars->n <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}
double birth_outcome(void * mypars)
{
	pars * my_pars = (pars *) mypars;
	++(my_pars->n) ; 
	return 0; 
}

void initial_conditions(void * mypars)
{
	pars * my_pars = (pars *) mypars;
	my_pars->n = No;
	my_pars->a = Ao;
	my_pars->time_index = 0;
	my_pars->hist_index = 0;
	my_pars->checkpts[0] = 0;
	my_pars->checkpts[1] = START_POLLUTING;
}




/** Execute all tasks that occur on fixed interval schedule 
  */
void fixed_interval_tasks(const double t, const void * mypars, void * myrecord)
{
	pars * my_pars = (pars *) mypars;
	int * my_record = (int *) myrecord;
	double * checkpts = my_pars->checkpts;
	

	/* Sample (print) system state at regular intervals */
	if(t > checkpts[0]){
		fprintf(stderr, "%lf %d\n", t, my_pars->n);
		checkpts[0] += SAMPLE_FREQ;
		my_record[my_pars->time_index] = my_pars->n;
		++my_pars->time_index;
	}

	/* Slowly change the bifurcation parameter */
	if(t > checkpts[1]){
		checkpts[1] += 1;  
		my_pars->a += 1;	
	}

}

void * pars_alloc()
{
	/** Allocates parameters structure.  
	 * Values that are reset each ensemble are set in 
	 * initial_conditions() function instead.  
	 */
	pars * my_pars = (pars *) malloc(sizeof(pars));
	my_pars->n = No;
	my_pars->K = 1000;
	my_pars->e = 0.5;
	my_pars->h = 200.0;
	return my_pars;
}

void popdyn(int * my_record)
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
	

	/** Parameters acutally allocated by par_alloc in gillespie call */
	void * my_pars = NULL; 
	size_t ensembles = ENSEMBLES;
	gillespie(rate_fn, outcome, N_EVENT_TYPES, my_pars, my_record, MAX_TIME, ensembles);
//	euler(my_pars, MAX_TIME, stderr);
//	gslode(my_pars, MAX_TIME, theory);
	free(my_pars);
}

int main(void){
	int my_record[1151];
	popdyn(my_record);
	return 0;
}

