#include "mach.h"
#include "vec_op.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

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
	if (size < 2 || size % 2 != 0)
	{
		printf("The number of processes needs to be a power of two, and at least two.\n");
		MPI_Finalize();
		return 2;
	}

	int n = atoi(argv[1]);
	double time_start, sum_1, sum_2, part_sum_1, part_sum_2, pi, error;
	double x_1 = (1.0 / 5);
	double x_2 = (1.0 / 239);
	double vec_1[n];
	double vec_2[n];
	int tag = 100;
	int d = n / (size - 1);
	double part_vec_1[n - (size - 2) * d];
	double part_vec_2[n - (size - 2) * d];

	if (rank == 0)
	{
		time_start = MPI_Wtime();
		compute_vec(vec_1, n, x_1);
		compute_vec(vec_2, n, x_2);
		for (int i = 1; i < size - 1; i++)
		{
			slice_vec(vec_1, part_vec_1, (i - 1) * d, i * d - 1);	
			slice_vec(vec_2, part_vec_2, (i - 1) * d, i * d - 1);
			MPI_Send(&part_vec_1, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);			
			MPI_Send(&part_vec_2, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);	
			MPI_Recv(&part_sum_1, 1, MPI_DOUBLE, i, tag + 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&part_sum_2, 1, MPI_DOUBLE, i, tag + 1, MPI_COMM_WORLD, &status);
			sum_1 = sum_1 + part_sum_1;
			sum_2 = sum_2 + part_sum_2;
			printf("sum: %f, %f\n", sum_1, sum_2);
		}
		slice_vec(vec_1, part_vec_1, (size - 2) * d, n - 1);
		slice_vec(vec_2, part_vec_2, (size - 2) * d, n - 1);
		MPI_Send(&part_vec_1, n, MPI_DOUBLE, size - 1, tag, MPI_COMM_WORLD);	
		MPI_Send(&part_vec_2, n, MPI_DOUBLE, size - 1, tag, MPI_COMM_WORLD);	
		MPI_Recv(&part_sum_1, 1, MPI_DOUBLE, size - 1, tag + 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&part_sum_2, 1, MPI_DOUBLE, size - 1, tag + 1, MPI_COMM_WORLD, &status);
		sum_1 = sum_1 + part_sum_1;
		sum_2 = sum_2 + part_sum_2;

		pi = 4 * (4 * sum_1 - sum_2);
		printf("Process %d pi: %f\n", rank, pi);
		printf("Walltime: %e (n = %d, np = %d)\n", MPI_Wtime() - time_start, n, size);
		error = fabs(M_PI - pi);
		printf("Error: %f (n = %d, np = %d)\n", error, n, size);
	}
	else
	{
		MPI_Recv(&part_vec_1, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&part_vec_2, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);	
		part_sum_1 = sum_vec(part_vec_1, n - (size - 2) * d);	
		part_sum_2 = sum_vec(part_vec_2, n - (size - 2) * d);	
		MPI_Send(&part_sum_1, 1, MPI_DOUBLE, 0, tag + 1, MPI_COMM_WORLD);		
		MPI_Send(&part_sum_2, 1, MPI_DOUBLE, 0, tag + 1, MPI_COMM_WORLD);	
	}	
	
	MPI_Finalize();
	return 0;
}
