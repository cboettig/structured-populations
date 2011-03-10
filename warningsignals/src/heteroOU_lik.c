/* **************************************** *
 *		ODE integration functions			*
 *	Provides GSL ode integrator and			*
 *	a simple euler integrator				*
 * **************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#define DIM 2


/* Define the ode system */
/* pars given as (alpha, m, theta, sigma) */
int func(double t, const double y[], double f[], void * mypars)
{
	double * pars = (double *) pars;
	double R = pars[0] + t*pars[1];
	f[0] = R*(pars[2] - y[1]); 
	f[1] = 2*R*y[1] + gsl_pow_2(pars[3]);
	return GSL_SUCCESS;
}


/* The goal is to avoid looping over the ode solver allocation */
/* Could write this a function that could be passed to an optimizer, 
 * should profile that to find out if it's worth passing the ode solver 
 * allocation along or just initializing each time.  */ 

void heteroOU(double * loglik, double * mypars,  double * X, double * times, int * N)
{

	/* Allocate space for mean and variance */
	double * Ex = (double *) malloc(*N * sizeof(double));
	double * Vx = (double *) malloc(*N * sizeof(double));


	/* Create our ODE system, we ignore the Jacobian */
	gsl_odeiv_system sys = {func, NULL, DIM, mypars};
	/* Define method as Embedded Runge-Kutta or RK Prince-Dormand (8,9) method */
	const gsl_odeiv_step_type * T
	 = gsl_odeiv_step_rk4;
//	 = gsl_odeiv_step_rk8pd;
	/* allocate stepper for our method of correct dimension*/
	gsl_odeiv_step * s
	 = gsl_odeiv_step_alloc (T, DIM);
	/* control will maintain absolute and relative error */
	gsl_odeiv_control * c
	 = gsl_odeiv_control_y_new (1e-6, 0.0);
	/* allocate the evolver */
	gsl_odeiv_evolve * e
	 = gsl_odeiv_evolve_alloc (DIM);

	/* Initial step size, will be modified as needed by adaptive alogorithm */
	double h = 1e-6;

	int i;
	for(i=0; i< *N; i++){
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




