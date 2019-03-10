#include "vec_op.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h> 

void print_vec(double *vec, int n)
{
	printf("Vector: ");
	for(int i = 0; i < n; i++)
	{
		printf("%f ", vec[i]);
	}
	printf("\n");
}

double sum_vec_using_openmp(double *vec, int len)
{
	double sum = 0;
	#pragma omp parallel for reduction(+:sum)
	for(size_t i = 0; i < len; i++)
	{
		sum = sum + vec[i];
	} 
	return sum;
}

void slice_vec(double *vec, double *part_vec, int start, int end)
{
	for(int i = 0; i <= end - start; i++)
	{
		part_vec[i] = vec[start + i];
	}
}




