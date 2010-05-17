#include "pars.h"

void * pars_alloc(int N, int K, double e, double a, double h, double start_polluting, double pollute_rate, double pollute_increment)
{

	pars * my_pars = (pars *) malloc(sizeof(pars));
	my_pars->n = N;
	my_pars->K = K;
	my_pars->e = e;
	my_pars->a = a;
	my_pars->h = h;
	my_pars->time_index = 0;
	my_pars->hist_index = 0;
	my_pars->checkpts[0] = 0;
	my_pars->checkpts[1] = my_pars->start_polluting;
	my_pars->start_polluting = start_polluting;
	my_pars->pollute_rate = pollute_rate;
	my_pars->pollute_increment = pollute_increment;
	return my_pars;
}
