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
	
	int n = atoi(argv[1]);
	int rank, size, i, tag;
	double sum_1, sum_2, pi;
	double x_1 = (1.0 / 5);
	double x_2 = (1.0 / 239);
	double vec_1[n];
	double vec_2[n];
	
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	tag = 100;
	if(rank == 0)
	{
		compute_vec(vec_1, n, x_1);
		compute_vec(vec_2, n, x_2);
		for(i = 1; i < size; i++)
		{
			MPI_Send(&vec_1, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
			MPI_Send(&vec_2, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);

		}
	}
	else
	{
		MPI_Recv(&vec_1, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&vec_2, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);	
		sum_1 = sum_vec(vec_1, n);
		sum_2 = sum_vec(vec_2, n);
        	pi = 4 * (4 * sum_1 - sum_2);
		printf("Process %d pi: %f\n", rank, pi);
	}
	MPI_Finalize();
	return 0;
}
