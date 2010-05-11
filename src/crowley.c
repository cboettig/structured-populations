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



/*		   0  1  2   3   4   5   6   7   8
 * Pars = {x, y, bx, by, dx, dy, cx, cy, K} */
void initial_conditions(void * ss)
{
	double * s = (double *) ss;
	s[0] = 50; 
	s[1] = 500;
	s[2] = 0.11;
	s[3] = 0.6;
	s[4] = 0.1;
	s[5] = 0.1;
	s[6] = 0.1;
	s[7] = 4.0;
	s[8] = 1000;
}


/** Execute all tasks that occur on fixed interval schedule 
 * Currently this function isn't very api like, as it is 
 * closely tied to the logical structure of the record data structure 
 * which itself could be abstracted more */
void fixed_interval_tasks(const double t, const void * mypars, void * myrecord)
{
	double * record = (double *) myrecord;
	if(t>record[0]){
		record[0] += 1;
		double * s = (double *) mypars;
		printf("%g %g %g\n", t, s[0], s[1]);
	}
}

/** record class for Crowley model */
typedef struct {
	double t_step;
	double * s1;
	double * s2;
} record;

record * = record_alloc(size_t N, double t_step)
{
	record * myrecord;
	myrecord->t_step = 10;
	myrecord->s1 = (double *) calloc(N, size(double));
	myrecord->s2 = (double *) calloc(N, size(double));
	return myrecord;
}

void record_free(record * myrecord)
{
	free(my_record->s1);
	free(my_record->s2);
	free(my_record);
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

	double my_record[3] = {0,0,0};
	double my_pars[9];
	initial_conditions(my_pars);

	

	const size_t n_event_types = 4;
	const size_t ensembles = 1;
	const size_t max_time = 500;
	gillespie(rate_fn, outcome, n_event_types, my_pars, my_record, max_time, ensembles);
//	euler(my_pars, MAX_TIME, stderr);
//	gslode(my_pars, MAX_TIME, theory);
}

int main(void)
{
	crowley();
	return 0;
}
