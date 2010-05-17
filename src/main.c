#include "warning_signals.h"
#define SAMPLE_TIME 1
#define ENSEMBLES 1
#define SAMPLE_FREQ 1.
#define MAX_TIME 500
#define START_POLLUTING 450
#define S 500 
void correlation(void);

int main(void){
	double  time[S], a[S], means[S], vars[S], skews[S], ar1[S], arN[S];
	double sample_freq = SAMPLE_FREQ, start_polluting = START_POLLUTING, pollute_rate = 1, pollute_increment = 1;
	int sample_time = SAMPLE_TIME, max_time = MAX_TIME, n_ensembles = ENSEMBLES;

	warning_signals(time, a, means, vars, skews, ar1, arN, &sample_time, &sample_freq, &max_time, &n_ensembles, &start_polluting, &pollute_rate, &pollute_increment);

	int i;
	for(i=0; i<S; i++){
		printf("%g %g %g\n", time[i], a[i], means[i] );
	}

	correlation();
	return 0;
}

