/**
 * @file crowley.c 
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
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


/*		   0  1  2   3   4   5   6   7   8
 * Pars = {x, y, bx, by, dx, dy, cx, cy, K} */

/** Specify all event functions.  Functions should return the rate as type double
 * and take the sole argument as a pointer to the pars structure.  */
double b1(void * vs)
{
	double * s = (double *) vs;
	      /* x *  bx *    K  -  x  -  y   +  cx  *   x  * y */
	return s[0] * s[2] * (s[8]-s[0]-s[1]) + s[6] * s[0] * s[1]  ; 
}
double b2(void * sv)
{ 
	double * s = (double *) vs;
	      /* y *  by *    K  -  x  -  y   +  cy  *   x  * y */
	return s[1] * s[3] * (s[8]-s[0]-s[1]) + s[7] * s[0] * s[1]  ; 
}
double d1(void * sv)
{ 
	double * s = (double *) vs;
	      /* x *  bx */
	return s[0] * s[4]  ; 
}

double d2(void * sv)
{ 
	double * s = (double *) vs;
	      /* x *  bx */
	return s[1] * s[5]  ; 
}


/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double b1_out(void * vs)
{ 
	double * s = (double *) vs;
	s[0] += 1;
}
double d1_out(void * vs)
{ 
	double * s = (double *) vs;
	s[0] -= 1;

	if(s[0] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}
double b2_out(void * vs)
{ 
	double * s = (double *) vs;
	s[1] += 1;
}
double d2_out(void * vs)
{ 
	double * s = (double *) vs;
	s[1] -= 1;

	if(s[0] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}



/*		   0  1  2   3   4   5   6   7   8
 * Pars = {x, y, bx, by, dx, dy, cx, cy, K} */
void initial_conditions(void * vs)
{
	double * s = (double *) vs;
	s[0] = 500; 
	s[1] = 500;
	s[2] = 0.2;
	s[3] = 0.6;
	s[4] = 0.1;
	s[5] = 0.1;
	s[6] = 0.1;
	s[7] = 0.1;
	s[8] = 1000;
}


/** Execute all tasks that occur on fixed interval schedule 
 * Currently this function isn't very api like, as it is 
 * closely tied to the logical structure of the record data structure 
 * which itself could be abstracted more */
void fixed_interval_tasks(const double t, const void * mypars, void * myrecord)
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

void crowley()
{
	/** Create a list of all event functions and their associated outcomes 
	 *  Order doesn't matter, but make sure outcomes are paired with the 
	 *  appropriate function! */
	event_fn rate_fn[4];	
	event_fn outcome[4];
	rate_fn[0] = &b1;
	outcome[0] = &b1_out;
	rate_fn[1] = &d1;
	outcome[1] = &d1_out;
	rate_fn[2] = &b2;
	outcome[2] = &b2_out;
	rate_fn[3] = &d2;
	outcome[3] = &d2_out;

	double * my_record = NULL;
	size_t ensembles = *n_ensembles;
	gillespie(rate_fn, outcome, 4, my_pars, my_record, *max_time, ensembles);
//	euler(my_pars, MAX_TIME, stderr);
//	gslode(my_pars, MAX_TIME, theory);
}


