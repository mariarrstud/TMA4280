#include "mach.h"
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
	double time_start, sum_1, sum_2, pi, error;
	double vec_1[n], vec_2[n];
	double x_1 = (1.0 / 5);
	double x_2 = (1.0 /239);

	time_start = omp_get_wtime();
	compute_vec(vec_1, n, x_1);
	compute_vec(vec_2, n, x_2);
	sum_1 = sum_vec_using_openmp(vec_1, n);
	sum_2 = sum_vec_using_openmp(vec_2, n);
	pi = 4 * (4 * sum_1 - sum_2);

	printf("pi: %f\n", pi);
	printf("Duration: %e (n = %d, nt = %d)\n", omp_get_wtime() - time_start, n, num_threads);
	
	error = fabs(M_PI - pi);
	printf("Error: %f (n = %d, nt = %d)\n", error, n, num_threads);
	
	return 0;
}
