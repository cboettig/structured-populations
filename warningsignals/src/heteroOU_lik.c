/* **************************************** *
 *		ODE integration functions			*
 *	Provides GSL ode integrator and			*
 *	a simple euler integrator				*
 * **************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#define DIM 2


/* Define the ode system */
/* pars given as (alpha, m, theta, sigma) */
int func(double t, const double y[], double f[], void *mypars)
{
	double * pars = (double *) mypars;
	double R = pars[0] + t*pars[1];
	f[0] = R*(pars[2] - y[1]); 
	f[1] = 2*R*y[1] + gsl_pow_2(pars[3]);
	return GSL_SUCCESS;
}


int
jac (double t, const double y[], double *dfdy,
  double dfdt[], void *mypars)
{
	double * pars = (double *) mypars;
	double sqrtR = sqrt( pars[0] + t*pars[1] );
	gsl_matrix_view dfdy_mat
	 = gsl_matrix_view_array (dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set (m, 0, 0, -2*sqrtR);
	gsl_matrix_set (m, 0, 1, 0.0);
	gsl_matrix_set (m, 1, 0, 0.0);
	gsl_matrix_set (m, 1, 1, -2*sqrtR);
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}


/* The goal is to avoid looping over the ode solver allocation */
/* Could write this a function that could be passed to an optimizer, 
 * should profile that to find out if it's worth passing the ode solver 
 * allocation along or just initializing each time.  */ 

void heteroOU(double *loglik, double *mypars,  double *X, double *times, int *N)
{

	/* Allocate space for mean and variance */
	double * Ex = (double *) malloc(*N * sizeof(double));
	double * Vx = (double *) malloc(*N * sizeof(double));


	/* Create our ODE system, we ignore the Jacobian */
	gsl_odeiv_system sys = {func, jac, DIM, mypars};
	/* Define method as Embedded Runge-Kutta Prince-Dormand (8,9) method */
	const gsl_odeiv_step_type * T
//	 = gsl_odeiv_step_rk8pd;
	 = gsl_odeiv_step_bsimp;
	/* allocate stepper for our method of correct dimension*/
	gsl_odeiv_step * s
	 = gsl_odeiv_step_alloc (T, DIM);
	/* control will maintain absolute and relative error */
	gsl_odeiv_control * c
	 = gsl_odeiv_control_y_new (1e-4, 1e-6);
	/* allocate the evolver */
	gsl_odeiv_evolve * e
	 = gsl_odeiv_evolve_alloc (DIM);

	/* Initial step size, will be modified as needed by adaptive alogorithm */
	double h = 1e-6;
	int i;

	
	for(i=0; i< (*N-1); i++){
		/* Range of time to simulate*/	
		double t = times[i], t1= times[i+1];
		/* initial conditions: start at X[i] with variance 0 */
		double y[DIM] = { X[i], 0.0 };

		/*dummy to make easy to switch to regular printing */
		double ti = t1; 
		while (t < ti)
		{
			int status = gsl_odeiv_evolve_apply (e, c, s,
												&sys,
												&t, t1,
												&h, y);
			if (status != GSL_SUCCESS)
			   break;
		}
		Ex[i] = y[0];
		Vx[i] = y[1];
	}

	*loglik = 0;
	/* first point X[0] observed at t[0], then prob X[1] given by Ex[0],Vx[0] */
	for(i=0; i < (*N-1); i++){
		*loglik += gsl_pow_2(Ex[i]-X[i+1]) / (2*Vx[i])  +.5*log(0.5*M_1_PI/Vx[i]);
	}


	free (Ex);
	free (Vx); 
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
}



int main(void)
{
	int i, N = 100;
	double pars[4] = {12., 0., 4., 8.};
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);
	double * X = (double *) malloc(N * sizeof(double));
	double * t = (double *) malloc(N * sizeof(double));
	for(i=0; i<N; i++){
		X[i] = gsl_ran_ugaussian(rng);
		t[i] = i;
	}
	double loglik = 0;

	heteroOU(&loglik, pars,  X, t, &N);
	printf("loglik = %e\n", loglik);
	free(X);
	free(t);
	gsl_rng_free(rng);
	return(0);
}
