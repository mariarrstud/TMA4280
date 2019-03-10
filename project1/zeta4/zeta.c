#include "zeta.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

void compute_vec(double *vec, int n)
{ 
	for(int i = 1; i <= n; i++)
	{
		vec[i - 1] = (1 / pow(i, 2));
	}
}
