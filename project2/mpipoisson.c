#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);
real solution(real x, real y);
void verification(real **u, size_t m, real *grid, real **error);
void inf_norm(real **error, size_t m);

void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
        return 1;
    }
    int n = atoi(argv[1]);
    if ((n & (n-1)) != 0) {
      printf("n must be a power-of-two\n");
      return 2;
    }
    int m = n - 1;
    real h = 1.0 / n;
    real *grid = mk_1D_array(n+1, false);
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }
    real *diag = mk_1D_array(m, false);
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }
    real **b = mk_2D_array(m, m, false);
    real **bt = mk_2D_array(m, m, false);
    int nn = 4 * n;
    real *z = mk_1D_array(nn, false);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);
        }
    }
    for (size_t i = 0; i < m; i++) {
        fst_(b[i], &n, z, &nn);
    }
    transpose(bt, b, m);
    for (size_t i = 0; i < m; i++) {
        fstinv_(bt[i], &n, z, &nn);
    }
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
        }
    }
    for (size_t i = 0; i < m; i++) {
        fst_(bt[i], &n, z, &nn);
    }
    transpose(b, bt, m);
    for (size_t i = 0; i < m; i++) {
        fstinv_(b[i], &n, z, &nn);
    }
    real **error = mk_2D_array(m, m, false);
    verification(b, m, grid, error);
    real norm = inf_norm(error);
    printf("Error: %e, h: %e\n", norm, h);
    return 0;
}

real rhs(real x, real y) {
    //return 1;
    verification_rhs = 5 * PI * PI * sin(PI * x) * sin(2 * PI * y);
    return verification_rhs;
}

void transpose(real **bt, real **b, size_t m)
{
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = b[j][i];
        }
    }
}

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    real **ret = (real **)malloc(n1 * sizeof(real *));
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

real solution(real x, real y)
{
	real sol = sin(PI * x) * sin(2 * PI * y);
	return sol;
}

void verification(real **u, real u_max, size_t m, real *grid, real **error)
{
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			error[i][j] = fabs(u[i][j] - solution(grid[i], grid[j]));
		}
	}
}

real inf_norm(real **error, size_t m)
{
	real norm = 0.0;
	for (size_t i = 0; i < m; i++)
	{
		real sum = 0.0;
		for (size_t j = 0; j < m; j++)
		{
			sum = sum + fabs(error[i][j]);
		}
		if (sum > norm)
		{
			norm = sum;
		}
	}
	return norm;
}
