#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>

double rhs(double x, double y)
{
	return 1.0;
}

void transpose(double *bt, double *b, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			bt[i][j] = b[j][i];
		}
	}
}

void fst(solution, *len_sol, support, *len_sup)
{

}

void fstinv(solution, *length1, z, *length2)
{

}

int poisson(int n)
{
	int m = n - 1;
	double h = 1.0 / n;
	int nn = 4 * n;
	int i, j;

	double grid[n + 1];
	for (i = 0; i < n + 1; i++)
	{
		grid[i] = i * h;
	}

	double diag[m]; //Storing eigenvalues of T
	for (i = 0; i < m; i++)
	{
		diag[i] = 2.0 * (1.0 - cos((i + 1) * M_PI / n));
	}
	double b[m]; //Storing U, G
	double bt[m]; //Storing ~U,~G
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			b[i][j] = h * h * rhs(grid[i + 1], grid[j + 1]);
		}
	}
	double z[nn];

	//1. ~G^T = S^-1 (S G)^T
	for (i = 0; i < m; i++)
	{
		fst(b[i], &n, z, &nn);
	}

	transpose(bt, b, m);
	for (i = 0; i < m; i++)
	{
		fstinv(bt[i], &n, z, &nn);
	}

	//2. ~U_ij = ~G_ij / (Lambda_i + Lambda_j)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
		}
	}

	//3. U = S^-1 (S ~U^T)^T
	for (i = 0; i < m; i++)
	{
		fst(bt[i], &n, z, &nn);
	}
	transpose(b, bt, m);
	for (i = 0; i < m; i++)
	{
		fstinv(b[i], &n, z, &nn);
	}

	double u_max = 0.0;
	for (i = 0; i  < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			u_max = umax > fabs(b[i][j]) ? u_max : fabs(b[i][j]);
		}
	}
	return 0;
}

