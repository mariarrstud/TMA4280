#include "zeta.h"
#include "vec_op.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("Invalid number of arguments\n");
	}
	int num_threads = atoi(argv[1]);
	omp_set_num_threads(num_threads);
	int n = atoi(argv[2]);
	double time_start, sum, pi, error;
	double vec[n];

	time_start = omp_get_wtime();
	compute_vec(vec, n);
	sum = sum_vec_using_openmp(vec, n);
	pi = sqrt(6 * sum);

	printf("pi: %f\n", pi);
	printf("Duration: %e (n = %d, nt = %d)\n", omp_get_wtime() - time_start, n, num_threads);
	
	error = fabs(M_PI - pi);
	printf("Error: %f (n = %d, nt = %d)\n", error, n, num_threads);
	
	return 0;
}
