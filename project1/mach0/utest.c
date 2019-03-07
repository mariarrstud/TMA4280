#include "mach.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv)
{
	double computed = compute_pi(3);
	double expected = 3.141621;
	if (abs(computed - expected) > 0)
	{
		printf("Testing failed\nExpected value: %f, Computed value: %f\n", expected, computed);
	}
	else
	{
		printf("Testing passed\n");
	}
	return 0;
}
