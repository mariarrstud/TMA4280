#include "mach.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

void compute_vec(double *vec, int n, double x)
{ 
	for(int i = 1; i <= n; i++)
	{
		vec[i - 1] = pow(- 1, i - 1) * pow(x, 2 * i - 1) * pow(2 * i - 1, - 1);
	}
}

