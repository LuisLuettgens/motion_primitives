#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "progressbar.hpp"

ProgressBar::ProgressBar(unsigned int maxval_) : maxval(maxval_) {
	memset(buffer, ' ', 76);
	size = 50;
	starttime = clock();
}

void ProgressBar::status(unsigned int current) {
	float percent = (float) current / (float) maxval * 100.f;
	set_buffer(percent);
	float time = time_elapsed();
	float torun = (100.f / percent * time) - time;
//	printf("%6.2f%% %s %7.2f seconds left\r", percent, buffer, torun);
	printf("%6d%% %s %7.1f sec left\r", (int)percent, buffer, torun);
	fflush(NULL);
}

ProgressBar::~ProgressBar() {
	set_buffer(100.0);
//	printf("%6.2f%% %s took %.2f seconds              \n", 100.0, buffer, time_elapsed());
	printf("%6d%% %s took %.1f sec              \n", 100, buffer, time_elapsed());
	fflush(NULL);
}


void ProgressBar::set_buffer(float percent) {
	buffer[0] = '[';
	unsigned int steps = (unsigned int) floor( percent / 100 * size);
	unsigned int pos = 1;
	for (; pos <= steps; pos++) {
		buffer[pos] = '-';
	}
	pos++;
	for (; pos < size; pos++) {
		buffer[pos] = ' ';
	}

	buffer[size] = ']';
	buffer[size +1] = '\0';
}

float ProgressBar::time_elapsed() {
	// save current time
	clock_t timer2 = clock();

	//printf("%d\n", CLOCKS_PER_SEC);
	float ret = (float) (timer2 - starttime) / (CLOCKS_PER_SEC);
	return ret;
}
