#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "zeta.h"

int main(int argc, char** argv)
{
	double computed = compute_pi(3);
	double expected = 2.857738;
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
