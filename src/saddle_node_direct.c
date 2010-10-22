#include "gillespie.h"
#include "gillespie_recording.h"
#define S_LENGTH 9
/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double sn_birth_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] += 1;
	return 0;
}
double sn_death_out(void * ss)
{ 
	double * s = (double *) ss;
	s[0] -= 1;

	if(s[0] <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}





/** Execute all tasks that occur on fixed interval schedule  */
void sn_fixed_interval(const double t, void * mypars, void * myrecord, int rep)
{
	record * my_record = (record *) myrecord;
	double * s = (double *) mypars;
	if (t > s[5]*my_record->t_step) 
	{
		my_record->s1[(int) s[5]+rep*my_record->N] = s[0]; 
//		printf("%g %g %g\n", t, s[0], s[2]);
		s[5] += 1; //increment sample counter

	/* increment the bifurcation parameter each sampletime after a certain wait */
	if(t > s[7]){ s[2] += s[6]*my_record->t_step; }
	}
}

/*		   0 1 2 3 4 5  6   7 8 
 * Pars = {n,e,a,K,h,i, Da, Dt, p} */
double sn_death(void * ss)
{
	const double * s = (double *) ss;
	/*		  n e    +  a      */
	return s[0]*s[1] + s[2]  ; 
}

double sn_birth(void * ss)
{ 
	const double * s = (double *) ss;
	      /*           e K n^2/ (n^2 +h^2)  */
		  int p = s[8];
	return gsl_pow_int(s[0], p)*s[1]*s[3] / ( gsl_pow_int(s[0],p) + gsl_pow_int(s[4],p) ); 
}


/** Must create a copy of parameter statespace which can be
 * modified in the loop.   */
void * sn_reset(const void * inits)
{
	const double * ss = (const double *) inits;
	double * s = (double *) calloc(S_LENGTH, sizeof(double));
	int i;
	for(i=0;i< S_LENGTH ;i++) s[i] = ss[i];
	s[5] = 0; // sampling counter
	return s;
}

void saddle_node_direct(double* s1, double* inits, int* n_samples, int* reps, double* maxtime)
{

	/** Create a list of all event functions and their associated outcomes 
	 *  Order doesn't matter, but make sure outcomes are paired with the 
	 *  appropriate function! */
	const size_t n_event_types = 2;
	event_fn rate_fn[n_event_types];	
	event_fn outcome[n_event_types];
	rate_fn[0] = &sn_birth;
	outcome[0] = &sn_birth_out;
	rate_fn[1] = &sn_death;
	outcome[1] = &sn_death_out;

	RESET reset_fn = &sn_reset;
	FIXED fixed_interval_fn = &sn_fixed_interval;	

	record * my_record = record_alloc(*n_samples, *reps, *maxtime);
	gillespie(rate_fn, outcome, n_event_types, inits, my_record, *maxtime, (size_t) *reps, reset_fn, fixed_interval_fn);

	int i;
	for(i = 0; i< (*n_samples)*(*reps); i++)
	{ 
		s1[i] = my_record->s1[i]; 
	}
	record_free(my_record);


/*	for(i=0; i< *n_samples; i++)
	{
		printf("%g\n",
			gsl_stats_mean( &(s1[i]),*n_samples, *reps)
			,gsl_stats_variance( &(s1[i]),*n_samples, *reps)
			);
	}
*/
}



int sn(void)
{
	const int n_samples = 100;
	const int replicates = 1;
	double s1[n_samples*replicates];


	//			0  1  2  3  4  5  6  7  8 
	//	Pars = {n, e, a, K, h, i, Da,Dt,p  } 
	double c_inits[9] = {572, .5, 160, 1000, 200, 0, 1, 100, 2.};

	double maxtime = 150; 
	int n = n_samples;
	int reps = replicates;

	saddle_node_direct(s1, c_inits, &n, &reps, &maxtime);

	printf("mean: mean = %g stdev = %g\n",
			gsl_stats_mean(s1, 1, n_samples), 
			sqrt(gsl_stats_variance(s1, 1, n_samples))
	); 
	return 0;
}
