#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

#include "vec_op.h"
#include "poisson.h"

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		printf("Invalid numbr of arguments\n");
		return 1;
	}

	int n = atoi(argv[1]);
	if ((n & (n - 1)) != 0)
	{
		printf("n must be a power of two\n");
		return 2;
	}
	
	poisson(n);	
	return 0;
}










