#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <time.h>
#include <omp.h>
#define _RNDSEED 1 /* 1 for true */

typedef double (* event_fn)(void * my_pars);
void * pars_cpy(void * in);

void fixed_interval_tasks(const double t, const void * my_pars, void * my_record);
void initial_conditions(void * mypars);

void gillespie( const event_fn * rate_fn, 
				const event_fn * outcome, 
				const size_t n_event_types, 
				void * my_pars,
				void * my_record,
				const size_t max_time, 
				const size_t ensembles 
				);
