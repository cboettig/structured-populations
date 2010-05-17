#include "warning_signals.h"
#define SAMPLE_TIME 50
#define ENSEMBLES 100
#define SAMPLE_FREQ 1.
#define MAX_TIME 500
#define START_POLLUTING 600

void correlation(void);

int main(void){
	double  time[1151], a[1151], means[1151], vars[1151], skews[1151], ar1[1151], arN[1151];
	double sample_freq = SAMPLE_FREQ, start_polluting = START_POLLUTING, pollute_rate = 1, pollute_increment = 1;
	int sample_time = SAMPLE_TIME, max_time = MAX_TIME, n_ensembles = ENSEMBLES;
	warning_signals(time, a, means, vars, skews, ar1, arN, &sample_time, &sample_freq, &max_time, &n_ensembles, &start_polluting, &pollute_rate, &pollute_increment);


	correlation();
	return 0;
}

