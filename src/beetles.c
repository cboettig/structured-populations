/**
 * @file beetles.c 
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 */

#include "gillespie.h"

typedef struct {
	int i;
	double t_step;
	double * s1;
	double * s2;
	double * s3;
	double * s4;
} record;

record* b_record_alloc(size_t N, double maxtime)
{
	record* myrecord = (record*) malloc(sizeof(record));
	myrecord->t_step = maxtime/N;
	myrecord->i = 0;
	myrecord->s1 = (double*) calloc(N+1, sizeof(double));
	myrecord->s2 = (double*) calloc(N+1, sizeof(double));
	myrecord->s3 = (double*) calloc(N+1, sizeof(double));
	myrecord->s4 = (double*) calloc(N+1, sizeof(double));
return myrecord;
}

void b_record_free(record * myrecord)
{
	free(myrecord->s1);
	free(myrecord->s2);
	free(myrecord->s3);
	free(myrecord->s4);
	free(myrecord);
}

/** Specify all event functions.  Functions should return 
 * the rate as type double and take the sole argument as 
 * a pointer to the pars structure.   */
double bE(void * ss)
{
	double * s = (double *) ss;
	return s[4]*s[3]; // b*A

}
double dE(void * ss)
{ 
	double * s = (double *) ss;
	return s[0]*(s[5]+s[3]*s[14] +s[1]*s[12]); // E*(ue+A*cae+L*cle)
}

double dL(void * ss)
{ 
	double * s = (double *) ss;
	return s[1]*s[6]; // L*ul 
}
double dP(void * ss)
{ 
	double * s = (double *) ss;
	return s[2]*(s[7] + s[3]*s[13]); // P*(up+A*cap)
}

double dA(void * ss)
{ 
	double * s = (double *) ss;
	return s[3]*s[8]; //A*ua
}

/** Maturation (two-transition event) */
double mE(void * ss)
{ 
	double * s = (double *) ss;
	return s[0]*s[9]; // E*ae
}

double mL(void * ss)
{ 
	double * s = (double *) ss;
	return s[1]*s[10]; //L*al
}

double mP(void * ss)
{ 
	double * s = (double *) ss;
	return s[2]*s[11]; // P*ap
}




/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double bE_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] += 1;
//	printf("birth, %g %g\n", s[3], bE(ss)); //error checking
	return 0;
}
double dE_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] -= 1;
	return 0 ; 
}
double dL_out(void * ss)
{ 
	double * s = (double *) ss;
	s[1] -= 1;
	return 0 ; 
}
double dP_out(void * ss)
{ 
	double * s = (double *) ss;
	s[2] -= 1;
	return 0 ; 
}
double dA_out(void * ss)
{ 
	double * s = (double *) ss;
	s[3] -= 1;
	return 0 ; 
}
double mE_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] -= 1;
	s[1] += 1;
	return 0 ; 
}
double mL_out(void * ss)
{ 
	double * s = (double *) ss;
	s[1] -= 1;
	s[2] += 1;
	return 0 ; 
}
double mP_out(void * ss)
{ 
	double * s = (double *) ss;
	s[2] -= 1;
	s[3] += 1;
	return 0 ; 
}





/** Must create a copy of parameter statespace which can be
 * modified in the loop.  Called from gillespie with no 
 * way to free this!  will create memory leak. must either make 
 * a free function or make this not allocate and instead pass two 
 * copies, one of which must become private */
void * beetles_reset(const void * inits)
{
	const double * ss = (const double *) inits;
	double * s = (double *) calloc(15, sizeof(double));
	int i;
	for(i=0;i<15;i++) s[i] = ss[i];
	return s;
}

/** Execute all tasks that occur on fixed interval schedule 
 * Currently this function isn't very api like, as it is 
 * closely tied to the logical structure of the record data structure 
 * which itself could be abstracted more */
void beetles_fixed_interval(const double t, const void * mypars, void * myrecord)
{
	record * my_record = (record *) myrecord;
	if (t > my_record->i * my_record->t_step) 
	{
		double * s = (double *) mypars;
		my_record->s1[my_record->i] = s[0]; 
		my_record->s2[my_record->i] = s[1]; 
		my_record->s3[my_record->i] = s[2]; 
		my_record->s4[my_record->i] = s[3]; 
		printf("%g %g %g %g %g\n", t, s[0], s[1], s[2], s[3]);
		++my_record->i;
	}
}

void beetles(double* s1, double* s2, double* s3, double* s4, double* inits, int* n_samples, double* maxtime)
{
	/** Create a list of all event functions and their associated outcomes 
	 *  Order doesn't matter, but make sure outcomes are paired with the 
	 *  appropriate function! */
	const size_t n_event_types = 8;
	event_fn rate_fn[n_event_types];	
	event_fn outcome[n_event_types];
	rate_fn[0] = &bE;
	outcome[0] = &bE_out;
	rate_fn[1] = &dE;
	outcome[1] = &dE_out;
	rate_fn[2] = &dL;
	outcome[2] = &dL_out;
	rate_fn[3] = &dP;
	outcome[3] = &dP_out;
	rate_fn[4] = &dA;
	outcome[4] = &dA_out;
	rate_fn[5] = &mE;
	outcome[5] = &mE_out;
	rate_fn[6] = &mL;
	outcome[6] = &mL_out;
	rate_fn[7] = &mP;
	outcome[7] = &mP_out;

	RESET reset_fn = &beetles_reset;
	FIXED fixed_interval_fn = &beetles_fixed_interval;

	const size_t ensembles = 1;
	
	record * my_record = b_record_alloc(*n_samples, *maxtime);

	gillespie(	rate_fn, outcome, n_event_types, 
				inits, my_record, *maxtime, 
				ensembles, reset_fn, fixed_interval_fn);


	int i;
	for(i = 0; i< *n_samples; i++)
	{ 
		s1[i] = my_record->s1[i]; 
		s2[i] = my_record->s2[i];
		s3[i] = my_record->s3[i]; 
		s4[i] = my_record->s4[i];
	}
	b_record_free(my_record);
}


int beetle(void)
{
	const int n_samples = 500;
	double s1[n_samples];
	double s2[n_samples];
	double s3[n_samples];
	double s4[n_samples];
//	       0  1  2  3  4   5   6   7   8   9  10   11  12   13  14  
// Pars = {E, L, P, A, b, ue, ul, up, ua, ae, al, ap, cle, cap, cae} 
	double inits[15] = {100, 0, 0, 0, 
						5., 
						0, 0.001, 0, 0.003,
						1/3.8, 1/(20.2-3.8), 1/(25.5-20.2), 
						0.01, 0.004, 0.01};

	int n = n_samples;
	double maxtime = 5000; 

	beetles(s1, s2, s3, s4, inits, &n, &maxtime);

	printf("m1 = %g sd1 = %g nm2 = %g sd2 = %g m3 = %g sd3 = %g m4 = %g sd4 = %g\n",
			gsl_stats_mean(s1, 1, n_samples), 
			sqrt(gsl_stats_variance(s1, 1, n_samples)),
			gsl_stats_mean(s2, 1, n_samples), 
			sqrt(gsl_stats_variance(s2, 1, n_samples)),
			gsl_stats_mean(s3, 1, n_samples), 
			sqrt(gsl_stats_variance(s3, 1, n_samples)),
			gsl_stats_mean(s4, 1, n_samples), 
			sqrt(gsl_stats_variance(s4, 1, n_samples))
	); 
	return 0;
}
