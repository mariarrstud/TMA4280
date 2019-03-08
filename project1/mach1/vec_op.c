#include "vec_op.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 

void print_vec(double *vec, int n)
{
	printf("Vector: ");
	for(int i = 0; i < n; i++)
	{
		printf("%f ", vec[i]);
	}
}

double sum_vec(double *vec, int len)
{
	double sum = 0;
	for(int i = 0; i < len; i++)
	{
		sum = sum + vec[i];
	}
	return sum;
}
