#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

double compute_vec(int n)
{
	double vec[n + 2]; 
	for(int i = 1; i <= n; i++)
	{
		vec[i] = (1 / pow(i, 2));
	}
	return vec;
}

int main(int argc, char **argv)
{
	int n = atoi(argv[1]);
	int rank, size;
	double sum, pi, vec[n + 2];
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank == 0)
	{
		vec = compute_vec(n);
	}
	else
	{
		MPI_Reduce(&vec, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	pi = sqrt(6 * sum);
	printf("pi%f\n", pi);	
	MPI_Finalize();
	return 0;
}
