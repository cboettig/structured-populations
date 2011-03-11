#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_siman.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200     	/* how many points do we try before stepping */
#define ITERS_FIXED_T 200	/* how many iterations for each T? */
#define STEP_SIZE .002 		/* max step size in random walk */
#define K 1.0				/* Boltzmann constant */
#define T_INITIAL 0.05		/* initial temperature */
#define MU_T 1.004		    /* damping factor for temperature */
#define T_MIN .008

/* setup parameters for multimin method*/
#define INIT_STEP 0.002
#define MAX_ITER 5000
#define ERR_TOL 1e-6
#define PRINT 1
#include <gsl/gsl_multimin.h>


double optim_func (const gsl_vector *v, void *params);

double multimin(gsl_vector *x, void * params);
double siman(gsl_vector * x, void * params, gsl_rng * rng); 

