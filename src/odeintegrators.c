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
#include "warning_signals.h"
#define DIM 2
#define NPTS 50.0 // number of points to print out, best as double notation

#define No 572

double round(double);

double death(void * mypars);
double birth(void * mypars);

double alpha1(double phi, void * mypars)
{
	pars * my_pars = (pars *) mypars;
	my_pars->n = (int) round(phi);
	return birth(mypars)-death(mypars); 
}
double alpha2(double phi, void * mypars)
{
	pars * my_pars = (pars *) mypars;
	my_pars->n = (int) round(phi);
	return birth(mypars) + death(mypars); 
}
double alpha1_prime(double phi, void * mypars)
{
	pars * my_pars = (pars *) mypars;
	my_pars->n = (int) round(phi);
	double v = (gsl_pow_2(my_pars->n) + gsl_pow_2(my_pars->h) );
	double vp= 2*my_pars->n;
	double u = my_pars->e * my_pars->K * gsl_pow_2( my_pars->n ); 
	double up = 2*my_pars->e * my_pars->K *  my_pars->n;
	return (up*v-u*vp)/gsl_pow_2(v) - my_pars->e; 
}

double alpha1_doubleprime(double phi, void * mypars)
{
	pars * my_pars = (pars *) mypars;
	my_pars->n = (int) round(phi);
	double v = (gsl_pow_2(my_pars->n) + gsl_pow_2(my_pars->h) );
	double h = my_pars->h;
	double n = my_pars->n;
	double e = my_pars->e;
	double K = my_pars->K;
	return 2*e*h*h*K*(3*n*n-h*h)/gsl_pow_3(v); 
}



/* Define the ode system */
int func(double t, const double y[], double f[], void * mypars)
{
	f[0] = alpha1(y[0], mypars); // + y[1]*alpha1_doubleprime(y[0], mypars)/2;
	f[1] = 2*y[1]*alpha1_prime(y[0], mypars) + alpha2(y[0], mypars);
	return GSL_SUCCESS;
}


/* Many good algorithms don't need the Jacobian, so we'll just pass a null pointer */
void gslode(void * mypars, double max_time, FILE *theory)
{
	/* Create our ODE system */
	gsl_odeiv_system sys = {func, NULL, DIM, mypars};

	/* Range of time to simulate*/	
	double t = 0.0, t1 = max_time;
	/* Initial step size, will be modified as needed by adaptive alogorithm */
	double h = 1e-6;
	/* initial conditions, No */
	double y[DIM] = { No, 0.0 };


	/* Define method as Embedded Runge-Kutta Prince-Dormand (8,9) method */
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

	/*dummy to make easy to switch to regular printing */
	double ti = t1; 
	int i;

	/* Uncomment the outer for loop to print *
	 * Npts sample pts at regular intervals   */
	for (i = 1; i <= NPTS; i++){   ti = i * t1 / NPTS;
		while (t < ti)
		{
			int status = gsl_odeiv_evolve_apply (e, c, s,
												&sys,
												&t, t1,
												&h, y);
			if (status != GSL_SUCCESS)
			   break;
//			fprintf(theory,"%.6e %.6e %.6e\n",t,y[0],y[1]); //adaptive mesh
		}
		fprintf(theory,"%g %g %g\n",t,y[0],y[1]);  }

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);
}


void euler(void *mypars, double max_time, FILE *theory){
 	
	double dt = .01;
	int sample; 

	double y[DIM] = { No, 0.0 };
	double f[DIM] = {0};
	double t=0.0;
	

	for(sample = 0; sample < max_time/dt; sample++){

		if( sample % (int) (1/dt) == 0)
		{

			fprintf
			(
				theory, "%g %g %g\n",
				dt*sample, 
				y[0], 
				y[1]
			);
		}
		int status = func(t, y, f, mypars);
		if(status != GSL_SUCCESS)
			break;
		y[0] +=  dt*f[0];
		y[1] +=  dt*f[1];
	}

}
