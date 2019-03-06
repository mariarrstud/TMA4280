#include "mach.h"

#include <stdio.h>
#include <math.h> 

double compute_pi(int n)
{
	double x_1 = (1.0 / 5);
	double x_2 = (1.0 / 239);
	double S_1 = 0;	
	double S_2 = 0;
	for(int i = 1; i <=  n; i = i + 1)
	{
		S_1 = S_1 + pow(- 1, i - 1) * pow(x_1, 2 * i - 1) * pow(2 * i - 1, - 1);
		S_2 = S_2 + pow(- 1, i - 1) * pow(x_2, 2 * i - 1) * pow(2 * i - 1, - 1);
	}
	double pi = 4 * (4 * S_1 - S_2);
	return pi;
}

