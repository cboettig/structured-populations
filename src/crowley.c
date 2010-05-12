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
#include "gillespie.h"

/** record class for Crowley model */
typedef struct {
	int i;
	double t_step;
	double * s1;
	double * s2;
} record;

record* record_alloc(size_t N, double maxtime)
{
	record* myrecord = (record*) malloc(sizeof(record));
	myrecord->t_step = maxtime/N;
	printf("%g\n", maxtime/N);
	myrecord->i = 0;
	myrecord->s1 = (double*) calloc(N+1, sizeof(double));
	myrecord->s2 = (double*) calloc(N+1, sizeof(double));
	return myrecord;
}

void record_free(record * myrecord)
{
	free(myrecord->s1);
	free(myrecord->s2);
	free(myrecord);
}



/*		   0  1  2   3   4   5   6   7   8
 * Pars = {x, y, bx, by, dx, dy, cx, cy, K} */

/** Specify all event functions.  Functions should return the rate as type double
 * and take the sole argument as a pointer to the pars structure.  */
double b1(void * ss)
{
	double * s = (double *) ss;
	      /* x *  bx *    K  -  x  -  y   +  cx  *   x  * y */
	return s[0] * s[2] *  (s[8]-s[0]-s[1])/s[8] + s[6] * s[0] * s[1]/s[8]  ; 
}
double b2(void * ss)
{ 
	double * s = (double *) ss;
	      /* y *  by *    K  -  x  -  y */
	return s[1] * s[3] * (s[8]-s[0]-s[1])/s[8]  ; 
}
double d1(void * ss)
{ 
	double * s = (double *) ss;
	      /* x *  bx */
	return s[0] * s[4]  ; 
}

double d2(void * ss)
{ 
	double * s = (double *) ss;
	      /* x *  bx   +   cy *  x   *  y  */
	return s[1] * s[5] + s[7] * s[0] * s[1]/s[8] ; 
}


/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double b1_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] += 1;
	return 0;
}
double d1_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] -= 1;

	if(s[0] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}
double b2_out(void * ss)
{ 
	double * s = (double *) ss;
	s[1] += 1;
	return 0;
}
double d2_out(void * ss)
{ 
	double * s = (double *) ss;
	s[1] -= 1;

	if(s[1] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}





/** Must create a copy of parameter statespace which can be
 * modified in the loop.  Will also take */
void * reset(void * inits)
{
	double * ss = (double *) inits;
	double * s = (double *) calloc(9, sizeof(double));
	int i;
	for(i=0;i<9;i++) s[i] = ss[i];
	return s;
}

/** Execute all tasks that occur on fixed interval schedule 
 * Currently this function isn't very api like, as it is 
 * closely tied to the logical structure of the record data structure 
 * which itself could be abstracted more */
void fixed_interval_tasks(const double t, const void * mypars, void * myrecord)
{
	record * my_record = (record *) myrecord;
	if (t > my_record->i * my_record->t_step) 
	{
		double * s = (double *) mypars;
		my_record->s1[my_record->i] = s[0]; 
		my_record->s2[my_record->i] = s[1]; 
		printf("%g %g %g\n", t, s[0], s[1]);
		++my_record->i;
	}
}

void crowley(double* s1, double* s2, double* inits, int* n_samples)
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

	
	const size_t n_event_types = 4;
	const size_t ensembles = 1;
	const size_t max_time = 5000;
	
	record * my_record = record_alloc(*n_samples, max_time);

	gillespie(rate_fn, outcome, n_event_types, inits, my_record, max_time, ensembles);

	int i;
	for(i = 0; i< *n_samples; i++)
	{ 
		s1[i] = my_record->s1[i]; 
		s2[i] = my_record->s2[i];
	}
	free(my_record);
//	euler(my_pars, MAX_TIME, stderr);
//	gslode(my_pars, MAX_TIME, theory);
}

int main(void)
{
	int n_samples = 100;
	double s1[100];
	double s2[100];
	/*			0  1  2   3   4   5   6   7   8
		Pars = {x, y, bx, by, dx, dy, cx, cy, K} */
	double inits[9];
	inits[0] = 500; 
	inits[1] = 4500;
	inits[2] = 0; //0.11;
	inits[3] = 0.6;
	inits[4] = 0; //0.1;
	inits[5] = 0.1;
	inits[6] = 0; // 0.1;
	inits[7] = 4.0;
	inits[8] = 10000;

	crowley(s1, s2, inits, &n_samples);
	printf("m1 = %g sd1 = %g\nm2 = %g sd2 = %g\n",
			gsl_stats_mean(s1, 1, n_samples), 
			sqrt(gsl_stats_variance(s1, 1, n_samples)),
			gsl_stats_mean(s2, 1, n_samples), 
			sqrt(gsl_stats_variance(s2, 1, n_samples))
	); 
	return 0;
}
