#include "zeta.h"
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
	double sum, pi;
	double vec[n];

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	tag = 100;
	if(rank == 0)
	{
		compute_vec(vec, n);
		for(i = 1; i < size; i++)
		{
			MPI_Send(&vec, n, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);

		}
	}
	else
	{
		MPI_Recv(&vec, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		sum = sum_vec(vec, n);
        	pi = sqrt(6 * sum);
		printf("Process %d pi: %f\n", rank, pi);	
	}	
	
	MPI_Finalize();
	return 0;
}
