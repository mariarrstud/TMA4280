#include "zeta.h"
#include "vec_op.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

int is_power_of_two(int number)
{
	while(number != 1)
	{
		if(number % 2 != 0)
		{
			return 0;
		}
		number = number / 2;
	}
	return 1;
}

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("Invalid number of arguments\n");
	}
	
	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	if (size < 2 || !is_power_of_two(size))
	{
		printf("The number of processes needs to be a power of two, and at least two.\n");
		MPI_Finalize();
		return 2;
	}
	
	int num_threads = atoi(argv[1]);
	omp_set_num_threads(num_threads);
	int n = atoi(argv[2]);
	double time_start, sum, part_sum, pi, error;
	double vec[n];
	int tag = 100;
	int d = n / (size - 1);
	double part_vec[n - (size - 2) * d];

	if (rank == 0)
	{
		time_start = MPI_Wtime();
		compute_vec(vec, n);
		for (int i = 1; i < size - 1; i++)
		{
			slice_vec(vec, part_vec, (i - 1) * d, i * d - 1);
			MPI_Send(&part_vec, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);	
			MPI_Recv(&part_sum, 1, MPI_DOUBLE, i, tag + 1, MPI_COMM_WORLD, &status);
			sum = sum + part_sum;
		}
		slice_vec(vec, part_vec, (size - 2) * d, n - 1);
		MPI_Send(&part_vec, n, MPI_DOUBLE, size - 1, tag, MPI_COMM_WORLD);	
		MPI_Recv(&part_sum, 1, MPI_DOUBLE, size - 1, tag + 1, MPI_COMM_WORLD, &status);
		sum = sum + part_sum;

		pi = sqrt(6 * sum);
		printf("Process %d pi: %f\n", rank, pi);
		printf("Duration: %e (n = %d, np = %d, nt = %d)\n", MPI_Wtime() - time_start, n, size, num_threads);
		error = fabs(M_PI - pi);
		printf("Error: %f (n = %d, np = %d, nt = %d)\n", error, n, size, num_threads);
	}
	else
	{
		MPI_Recv(&part_vec, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		part_sum = sum_vec_using_openmp(part_vec, n - (size - 2) * d);	
		MPI_Send(&part_sum, 1, MPI_DOUBLE, 0, tag + 1, MPI_COMM_WORLD);	
	}	
	
	MPI_Finalize();
	return 0;
}
