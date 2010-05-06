#include "kde.h"

double brd0(double x[], const int N)
{
	gsl_sort(x, 1, N);
	double hi = gsl_stats_sd(x, 1, N);
	double iqr =
		gsl_stats_quantile_from_sorted_data (x,1, N,0.75) - 
        gsl_stats_quantile_from_sorted_data (x,1, N,0.25);
	double lo = GSL_MIN(hi, iqr/1.34);
	double bw = 0.9 * lo * pow(N,-0.2);
	return(bw);
}

/* kernels for kernel density estimates */
double gauss_kernel(double x)
{ 
	return exp(-(gsl_pow_2(x)/2))/(M_SQRT2*sqrt(M_PI)); 
}

double kerneldensity(double *samples, double obs, size_t n)
{
	size_t i;
	double h = GSL_MAX(brd0(samples, n), 1e-6);
	double prob = 0;
	for(i=0; i < n; i++)
	{
		prob += gauss_kernel( (samples[i] - obs)/h)/(n*h);
	}
	return prob;
}


