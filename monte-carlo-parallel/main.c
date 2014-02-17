// Aimal Khan - 0946636 - SE4F03 Assignment 1 

// main.c

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

extern double fcn(double *x, int n);
extern double parallelMonteCarlo (double *a, double *b, int n, long int N, double (*fcn)(double *x, int n), int my_rank, int p, MPI_Comm com);
extern double checkResult(double integral, double *a, double *b, int n, long int N, double (*fcn)(double *x, int n));

int main(int argc, char *argv[])
{
	// MPI Variables
	int my_rank;
	int p;
	// Integral value to be returned
	double integral;

	// # of dimensions
	int n = atof(argv[1]);
	// # of integration points
	int N = atof(argv[argc - 1]);

	// Arrays
	double *a = (double *) malloc(n*sizeof(double));
	double *b = (double *) malloc(n*sizeof(double));

	// intake input parameters
	for (int i = 0; i < n; i++)
	{
		a[i] = atof(argv[i + 2]);
		b[i] = atof(argv[i + n + 2]);
	}

	// Initialize MPI
	MPI_Init(&argc, &argv);
	// Get process rank
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	// Get # of processes
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	integral = parallelMonteCarlo(a, b, n, N, &fcn, my_rank, p, MPI_COMM_WORLD);
	if (my_rank == 0)
	{
		checkResult(integral, a, b, n, N, &fcn);
		printf("Value of the integral is %.4e\n", integral);
	}

	// free arrays and exit MPI
	MPI_Finalize();
	free(a);
	free(b);
	return 0;
}