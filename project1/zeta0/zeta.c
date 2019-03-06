#include "zeta.h"

#include <stdio.h>
#include <math.h> 

double compute_pi(int n)
{
	double S = 0;
	for(int i = 1; i <=  n; i = i + 1)
	{
		S = S + (1 / pow(i, 2));
	}
	double pi = sqrt(6 * S);
	return pi;
}

