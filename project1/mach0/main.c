#include "mach.h"
#include "utest.h"
#include "vtest.h"

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
	int run_or_test = atoi(argv[1]);
	int n = atoi(argv[2]);
	if(run_or_test == 1)
	{
		double pi = compute_pi(n);
		printf("Approximation of pi for n = %d: %f\n", n, pi);
	}
	else if((run_or_test) == 2)
	{
		print_utest();
	}
	else
	{
		save_vtest();
	}
	return 0;
}
