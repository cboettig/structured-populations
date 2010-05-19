#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>

#define L 8000
#define K 1
#define R .09
#define DT 0.01
#define GAMMA 1
void acorr(double * data, double * corr, size_t N)
{
		int m, n;
		for(m=0; m<N; m++)
		{
			for(n=0; n< N-m-1; n++)
			{
				corr[m] +=  data[n+m]*data[n];
			}
				corr[m] = (double) corr[m]/N;
		}
}

void ar1_series(gsl_rng * rng, double * data, size_t N, double rho)
{
	int i;
	for(i=1; i<N; i++)
	{
		data[i] = rho*data[i-1] + gsl_ran_gaussian(rng, 1) * sqrt(1-gsl_pow_2(rho)) ;
	}
}

void ou1(gsl_rng * rng, double * data, size_t N, double alpha, double dt)
{
	int i;
	double theta = 0, sigma = 1;
	printf("# Theory: Var = %lf\n", log(.5*gsl_pow_2(sigma)/alpha));
	printf("# lambda = %lf\n", alpha);
	data[0] = 0;
	for(i=1; i<N; i++){
		data[i] = -alpha*(theta - data[i-1])*dt + sqrt(dt)*gsl_ran_gaussian(rng, sigma);
	}
}

void langevin(gsl_rng * rng, double * data, size_t N, double rho, double dt)
{
	int i;
	double KbT = 1, gamma = GAMMA, kappa = R;
	printf("# Theory: Var = %lf\n", log(KbT/kappa));
	printf("# lambda = %lf\n", kappa/gamma);

	data[0] = 0;
	for(i=1; i<N; i++){
		data[i] = data[i-1]*(1-kappa*dt/gamma) + dt*sqrt(2*KbT*gamma)*gsl_ran_gaussian(rng, 1)/gamma;
	}
}




void correlation(void)
{
	double dt = DT;

	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default);
	gsl_rng_set(rng, time(NULL));
	double data[L], corr[L], logcorr[L], time[L];
	int i;
	FILE * file = fopen("test.txt", "w");

	data[0] = 0; 	corr[0] = 0; time[0] = 0;

	for(i=1; i < L; i++)
	{
		corr[i] = 0;
		time[i] = i*dt;
	}


//	ar1_series(rng, data, L, R);
	//ou1(rng, data, L, R, dt);
	langevin(rng, data, L, R, dt);

	acorr(data, corr, L);
	for(i=0; i < L; i++)
	{
		if(corr[i] > 0){
			logcorr[i] = log(corr[i]/K);
		} else{  
			break;
		}
		fprintf(file, "%g %g %g\n", time[i], corr[i], data[i]);
	}
	i = i/6;


	double c0, c1, cov00, cov01, cov11, chisq;

	gsl_fit_linear (time, 1, logcorr, 1, i,
										&c0, &c1, &cov00, &cov01, &cov11,
										&chisq);

	printf ("# best fit: Y = %g + %g X\n", c0, c1);
	printf ("# covariance matrix:\n");
	printf ("# [ %g, %g\n#   %g, %g]\n",
				 cov00, cov01, cov01, cov11);
	printf ("# chisq = %g\n", chisq);
	printf ("# pts = %d\n", i );

	fclose(file);
}


