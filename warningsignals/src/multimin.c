#include "optimizers.h"

/*Function being optimized is optim_func */

double multimin(gsl_vector *x, void * params)
{
	int i;
	size_t iter = 0;
	int status;
	double size;
	size_t n_pars = x->size;

	/* Declare minimizer type and allocate it */
	gsl_multimin_fminimizer * s = 
		gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, n_pars);

	/* Set initial step sizes to a vector of INIT_STEP size) */
	gsl_vector *ss = gsl_vector_alloc (n_pars);
	gsl_vector_set_all (ss, INIT_STEP);

	/* Initialize method and iterate */
	gsl_multimin_function minimize_this;
	minimize_this.n = n_pars;  // dimension
	minimize_this.f = &optim_func;
	minimize_this.params = params;

	gsl_multimin_fminimizer_set (s, &minimize_this, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, ERR_TOL);



		if(PRINT){ // Printing turned off
			printf ("iter: %5d, llik = %3.3f, size = %.3f,\n par values = ",
				   iter,
				   -s->fval, 
				   size );
			for(i = 0; i < n_pars; i++)
				printf (" %1lf ", gsl_vector_get (s->x, i) );
			printf("\n");
		}


		if (status == GSL_SUCCESS)
		{
//			printf ("converged to minimum at\n");
		}
			
	} while (status == GSL_CONTINUE && iter < MAX_ITER);

	/* Print out status, whether or not it converged */
	size = gsl_multimin_fminimizer_size (s);
	if(PRINT){ // Printing turned off
		printf ("iter: %5d, llik = %3.3f, size = %.3f,\n par values = ",
			   iter,
			   -s->fval, 
			   size );
		for(i = 0; i < n_pars; i++)
			printf (" %1lf ", gsl_vector_get (s->x, i) );
		printf("\n");
	}
	double loglik = -s->fval;
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return(loglik);
}

