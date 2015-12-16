#include <omp.h>
#include <stdio.h>
#include <time.h>
#include "ttas_sleep.lock.c"
int main() {
	int i = 0, j = 0;
	long t0 = clock();
        omp_set_num_threads(4);
#pragma omp parallel for shared(i), private(j)
		for (j = 0; j < 200000000; j++) {		
			lock();
			i++;
			unlock();
		}
	t0 = clock() - t0;
	printf("Total count = %d in %ld ticks\n", i, t0);
}