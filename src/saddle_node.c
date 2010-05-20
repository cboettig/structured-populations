#include "warning_signals.h"


/** Specify all event functions.  Functions should return the rate as type double
 * and take the sole argument as a pointer to the pars structure.  */
double death(void * mypars)
{ 
	pars * my_pars = (pars *) mypars;
	return my_pars->n * my_pars->e + my_pars->a; 
}
double birth(void * mypars)
{ 
	pars * my_pars = (pars *) mypars;
	return my_pars->e * my_pars->K * gsl_pow_2( my_pars->n ) / 
			(gsl_pow_2(my_pars->n) + gsl_pow_2(my_pars->h) ); 
}



/** outcomes functions can return a flag of 1 to break the time simulation (i.e. if extinction occurs)
 *  They take an argument of type pars and update the state directly via this structure.  Otherwise
 *  they should return 0 for success.  */
double death_outcome(void * mypars)
{ 
	pars * my_pars = (pars *) mypars;
	--(my_pars->n);
	if(my_pars->n <= 0){
		fprintf(stderr, "extinction\n"); 
		return 1 ;
	}
	return 0 ; 
}
double birth_outcome(void * mypars)
{
	pars * my_pars = (pars *) mypars;
	++(my_pars->n) ; 
	return 0; 
}


