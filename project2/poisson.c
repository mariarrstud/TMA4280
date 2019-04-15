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
	int counts[size] = { [0 ... size] = rows_p};
	int displs[size] = { 0 };
	for (size_t i = 1; i < size; i++) {
		if (rem > 0) {
			displs[i] = displs[i - 1] + rows_p + 1;
			counts[i - 1] ++;
			rem --;
		}
		else {
			displs[i] = displs[i - 1] + rows_p;	
		}
	}
	
	printf("counts:\n");
	print_vec(counts, size);
	printf("displs:\n");
	print_vec(displs, size);
	
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
			b[i][j] += h * h * rhs(grid[i+1], grid[j+1]);
		}
	}
	//#pragma omp parallel for schedule(static) reduction(+: b, z)
	for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
		fst_(b[i], &n, z, &nn);
	}
	
	printf("Matrix b process %d :", rank);
	print_matrix(b, m);
	
	//Pack data into sendbuffer
	double sendbuf1[m * counts[1]];
	size_t ind_send1 = 0;
	for (size_t k = 0; k < size; k++) {
		for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
			for (size_t j = displs[k]; j < displs[k] + counts[k]; j++) {
				sendbuf1[ind_send1] = b[i][j];
				ind_send1 ++;
			}
		}
	}
	
	printf("Sendbuffer process %d :", rank);
	print_vec(sendbuf1, m  * counts[1]);
	
	double recvbuf1[m * counts[1]];
	//MPI_Alltoallv
	MPI_Alltoallv(&sendbuf1, counts, displs, MPI_Double, &recvbuf1, 
		      counts, displs, MPI_Double, MPI_Comm comm);
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
	
	printf("Recvbuffer process %d :", rank);
	print_vec(recvbuf1, m  * counts[1]);
	
	printf("Matrix bt process %d :", rank);
	print_matrix(bt, m);
	
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
	double sendbuf2[m * counts[1]];
	size_t ind_send2 = 0;
	for (size_t k = 0; k < size; k++) {
		for (size_t i = displs[rank]; i < displs[rank] + counts[rank]; i++) {
			for (size_t j = displs[k]; j < displs[k] + counts[k]; j++) {
				sendbuf2[ind_send2] = bt[i][j];
				ind_send2 ++;
			}
		}
	}
	double recvbuf2[m * counts[1]];
	//MPI_Alltoallv
	MPI_Alltoallv(sendbuf, counts, displs, MPI_Double, recvbuf, 
		      counts, displs, MPI_Double, MPI_Comm comm);
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
	
	double u_max = 0.0;
    	for (size_t i = 0; i < m; i++) {
        	for (size_t j = 0; j < m; j++) {
        		u_max = u_max > fabs(b[i][j]) ? u_max : fabs(b[i][j]);
        	}
    	}
	printf("u_max = %e\n", u_max);
	//printf("T%e: %e\n", size, duration);
	
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
	return 1;
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
