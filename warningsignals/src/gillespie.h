#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <time.h>
#include <omp.h>
#define _RNDSEED 1 /* 1 for true */

typedef double (* event_fn)(void * my_pars);
typedef void * (* RESET)(const void * inits);
typedef void (* FIXED)(const double t, void * my_pars, void * my_record, int rep);

void * reset(void * inits);

//void fixed_interval_tasks(const double t, const void * my_pars, void * my_record);



void gillespie( const event_fn * rate_fn, 
				const event_fn * outcome, 
				const size_t n_event_types, 
				void * my_pars,
				void * my_record,
				const size_t max_time, 
				const size_t ensembles,
				RESET reset_fn,
				FIXED fixed_interval_fn
				);
