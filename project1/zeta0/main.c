#include "zeta.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		printf("Invalid number of arguments.\n");
		return 1;
	}
	int n = atoi(argv[1]);
	double pi = compute_pi(n);
	printf("Approximation of pi for n = %d: %f\n", n, pi);
	return 0;
}
