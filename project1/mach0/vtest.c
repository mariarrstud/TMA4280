#include "mach.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void save_vtest()
{
	FILE *vtest = fopen("vtest.txt", "w");
	double error;
	if(vtest == NULL)
	{
		printf("Could not open file\n");
	}
	for(int k = 1; k <= 24; k = k + 1)
	{
		error = fabs(M_PI - compute_pi(pow(2, k)));
		fprintf(vtest, "Error for n = %d: %f\n", (int)pow(2, k), error);
	}
}

