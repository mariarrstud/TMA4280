#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, size_t m);
real rhs(real x, real y);
void print_vec(double *vec, int len);
void print_vec_int(int *vec, int len);
void print_matrix(double **mat, int len);

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
	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	if (size < 2 || (size & (size - 1)) != 0) {
		printf("The number of processes needs to be a power-of-two, and at least two\n");
		MPI_Finalize();
		return 3;
	}
	real h = 1.0 / n;
	int m = n - 1;
	if (size > m) {
		size = m;
	}
	int rows_p = m / size;
	int rem = m % size;
	int counts[size];
	counts[size - 1] = rows_p;
	int displs[size];
	displs[0] = 0;
	for (size_t i = 1; i < size; i++) {
		if (rem > 0) {
			displs[i] = displs[i - 1] + rows_p + 1;
			counts[i - 1] = rows_p + 1;
			rem --;
		}
		else {
			displs[i] = displs[i - 1] + rows_p;
			counts[i - 1] = rows_p;
		}
	}
	int mpi_counts[size];
	int mpi_displs[size];
	for (size_t i = 0; i < size; i++) {
		mpi_counts[i] = counts[i] * counts[rank];
		mpi_displs[i] = displs[i] * counts[rank];
	}
	
	real time_start = MPI_Wtime();
	
	real *grid = mk_1D_array(n + 1, false);
	//#pragma omp parallel for schedule(static) reduction(+: grid)
	for (size_t i = 0; i < n+1; i++) {
		grid[i] += i * h;
	}
	real *diag = mk_1D_array(m, false);
	//#pragma omp parallel for schedule(static) reduction(+: diag)
	for (size_t i = 0; i < m; i++) {
		diag[i] += 2.0 * (1.0 - cos((i+1) * PI / n));
	}
	real **b = mk_2D_array(m, m, false);
	real **bt = mk_2D_array(m, m, false);
	int nn = 4 * n;
	real *z = mk_1D_array(nn, false);
	//#pragma omp parallel for schedule(static) reduction(+: b)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i][j] = h * h * rhs(grid[i+1], grid[j+1]);
		}
	}
	
	//#pragma omp parallel for schedule(static) reduction(+: b, z)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		fst_(b[i], &n, z, &nn);
	}
	
	//Pack data into sendbuffer
	double sendbuf1[m * counts[rank]];
	size_t ind_send1 = 0;
	for (size_t k = 0; k < size; k++) {
		for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
			for (size_t j = displs[k]; j < displs[k] + counts[k]; j++) {
				sendbuf1[ind_send1] = b[i][j];
				ind_send1 ++;
			}
		}
	}
	
	double recvbuf1[m * counts[rank]];
	//MPI_Alltoallv
	MPI_Alltoallv(&sendbuf1, mpi_counts, mpi_displs, MPI_DOUBLE, &recvbuf1, 
		      mpi_counts, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
	//Unwrap data
	size_t ind_recv1 = 0;
	for (size_t k = 0; k < size; k++) {
		for (size_t j = displs[k]; j < displs[k] + counts[k]; j++) {
			for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
				bt[i][j] = recvbuf1[ind_recv1];
				ind_recv1 ++;
			}
		}
	}
	
	//#pragma omp parallel for schedule(static) reduction(+: bt, z)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		fstinv_(bt[i], &n, z, &nn);
	}
	//#pragma omp parallel for schedule(static) reduction(+: bt)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		for (size_t j = 0; j < m; j++) {
			bt[i][j] = bt[i][j] / (diag[i] + diag[j]);
		}
	}
	//#pragma omp parallel for schedule(static) reduction(+: bt, z)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		fst_(bt[i], &n, z, &nn);
	}
	
	//Pack data into sendbuffer
	double sendbuf2[m * counts[rank]];
	size_t ind_send2 = 0;
	for (size_t k = 0; k < size; k++) {
		for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
			for (size_t j = displs[k]; j < displs[k] + counts[k]; j++) {
				sendbuf2[ind_send2] = bt[i][j];
				ind_send2 ++;
			}
		}
	}
	double recvbuf2[m * counts[rank]];
	//MPI_Alltoallv
	MPI_Alltoallv(&sendbuf2, mpi_counts, mpi_displs, MPI_DOUBLE, &recvbuf2, 
		      mpi_counts, mpi_displs, MPI_DOUBLE, MPI_COMM_WORLD);
	//Unwrap data
	size_t ind_recv2 = 0;
	for (size_t k = 0; k < size; k++) {
		for (size_t j = displs[k]; j < displs[k] + counts[k]; j++) {
			for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
				b[i][j] = recvbuf2[ind_recv2];
				ind_recv2 ++;
			}
		}
	}
	
	//#pragma omp parallel for schedule(static) reduction(+: b, z)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		fstinv_(b[i], &n, z, &nn);
	}
	
	real duration = time_start - MPI_Wtime();
	
	double local_u_max = 0.0;
	double u_max;
    	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
        	for (size_t j = 0; j < m; j++) {
        		local_u_max = local_u_max > fabs(b[i][j]) ? local_u_max : fabs(b[i][j]);
        	}
    	}
	MPI_Reduce(&local_u_max, &u_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		printf("u_max = %e\n", u_max);
		//printf("T%e: %e\n", size, duration);
	}
	
	MPI_Finalize();
	return 0;
}	

void print_vec(double *vec, int len)
{
	for (int i = 0; i < len; i++)
	{
		printf("%e ", vec[i]);
	}
	printf("\n");
}

void print_vec_int(int *vec, int len)
{
	for (int i = 0; i < len; i++)
	{
		printf("%d ", vec[i]);
	}
	printf("\n");
}

void print_matrix(double **mat, int len)
{
	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < len; j++)
		{
			printf("%e ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

real rhs(real x, real y) {
	//return 1;
	return (5 * PI * PI * sin(PI * x) * sin(2 * PI *y));
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
