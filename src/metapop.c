/**
 * @file metapop.c 
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
#include "gillespie.h"
#include "gillespie_recording.h"

/** Execute all tasks that occur on fixed interval schedule  */
void m_fixed_interval(const double t, void * mypars, void * myrecord, int rep)
{
	record * my_record = (record *) myrecord;
	double * s = (double *) mypars;
	if (t > s[9]*my_record->t_step) 
	{
		my_record->s1[(int) s[9]+rep*my_record->N] = s[0]; 
		my_record->s2[(int) s[9]+rep*my_record->N] = s[1]; 
		//printf("%g %g %g\n", t, s[0], s[1]);
		s[9] += 1; //increment sample counter
	}
}


/*		   0  1  2   3   4   5   6   7   8
 * Pars = {x, y, bx, by, dx, dy, cx, cy, K} */

/** Specify all event functions.  Functions should return the rate as type double
 * and take the sole argument as a pointer to the pars structure.  */
double bx(void * ss)
{
	const double * s = (double *) ss;
	/*		    x*bx*(K  -  x  -  y)    */
	return s[0]*s[2]*(s[8]-s[0]-s[1])/s[8]  ; 
}
double by(void * ss)
{ 
	const double * s = (double *) ss;
	      /* y *  by *    K  -  x  -  y */
	return s[1]*s[3]*(s[8]-s[0]-s[1])/s[8]  ; 
}
double dx(void * ss)
{ 
	const double * s = (double *) ss;
	      /* x *  bx */
	return s[0]*s[4]  ; 
}

double dy(void * ss)
{ 
	const double * s = (double *) ss;
	      /* y *  by   */
	return s[1]*s[5] ; 
}

double m(void * ss)
{ 
	const double * s = (double *) ss;
	      /*   cy*x*y	  */
	return  s[7]*s[0]*s[1]/s[8] ; 
}


/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double bx_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] += 1;
	return 0;
}
double dx_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] -= 1;

	if(s[0] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}
double by_out(void * ss)
{ 
	double * s = (double *) ss;
	s[1] += 1;
	return 0;
}
double dy_out(void * ss)
{ 
	double * s = (double *) ss;
	s[1] -= 1;

	if(s[1] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}

double m_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] += 1;
	s[1] -= 1;
	return 0;
}



/** Must create a copy of parameter statespace which can be
 * modified in the loop.  Will also take */
void * m_reset(const void * inits)
{
	const double * ss = (const double *) inits;
	double * s = (double *) calloc(10, sizeof(double));
	int i;
	for(i=0;i<9;i++) s[i] = ss[i];
	s[9] = 0; // sampling counter
	return s;
}

void metapop(double* s1, double* s2, double* inits, int* n_samples, int* reps, double* maxtime)
{

	/** Create a list of all event functions and their associated outcomes 
	 *  Order doesn't matter, but make sure outcomes are paired with the 
	 *  appropriate function! */
	const size_t n_event_types = 5;
	event_fn rate_fn[n_event_types];	
	event_fn outcome[n_event_types];
	rate_fn[0] = &bx;
	outcome[0] = &bx_out;
	rate_fn[1] = &dx;
	outcome[1] = &dx_out;
	rate_fn[2] = &by;
	outcome[2] = &by_out;
	rate_fn[3] = &dy;
	outcome[3] = &dy_out;
	rate_fn[4] = &m;
	outcome[4] = &m_out;

	RESET reset_fn = &m_reset;
	FIXED fixed_interval_fn = &m_fixed_interval;
	
	record * my_record = record_alloc(*n_samples, *reps, *maxtime);
	gillespie(rate_fn, outcome, n_event_types, inits, my_record, *maxtime, (size_t) *reps, reset_fn, fixed_interval_fn);

	int i;
	for(i = 0; i< (*n_samples)*(*reps); i++)
	{ 
		s1[i] = my_record->s1[i]; 
		s2[i] = my_record->s2[i];
	}
	record_free(my_record);

/*
	for(i=0; i< *n_samples; i++)
	{
		printf("%g %g %g %g\n",
			gsl_stats_mean( &(s1[i]),*n_samples, *reps),
			gsl_stats_mean( &(s2[i]),*n_samples, *reps),
			gsl_stats_variance( &(s1[i]),*n_samples, *reps),
			gsl_stats_variance( &(s2[i]),*n_samples, *reps));
	}
*/
}



int meta(void)
{
	const int n_samples = 100;
	const int replicates = 20;
	double s1[n_samples*replicates];
	double s2[n_samples*replicates];

//	double mean[n_samples];
//	double var[n_samples];

	//			0  1  2   3   4   5   6   7   8
	//	Pars = {x, y, bx, by, dx, dy, cx, cy, K} 
	double c_inits[9] = {500, 500, 0.2, .6, .1, .1, .1, .1, 10000};

	double maxtime = 300; 
	int n = n_samples;
	int reps = replicates;

	metapop(s1, s2, c_inits, &n, &reps, &maxtime);

//	int i;
//	for(i = 0; i< n_samples * replicates; i++)
//		printf("%g, %g\n", s1[i], s2[i] );

	printf("metapop: m1 = %g sd1 = %g m2 = %g sd2 = %g\n",
			gsl_stats_mean(s1, 1, n_samples), 
			sqrt(gsl_stats_variance(s1, 1, n_samples)),
			gsl_stats_mean(s2, 1, n_samples), 
			sqrt(gsl_stats_variance(s2, 1, n_samples))
	); 
	return 0;
}
