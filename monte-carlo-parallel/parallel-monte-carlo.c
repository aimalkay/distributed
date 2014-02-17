// Aimal Khan - 0946636 - SE4F03 Assignment 1 

// parallel-monte-carlo.c

// returns approximation for monte-carlo on process 0
// a is a pointer to n doubles, where a[i] stores a_i+1
// b is a pointer to n doubles, where a[i] stores b_i+1
// n is dimension n
// N is # of random points to be used in the integration
// fcn is a pointer to function returning double with arguments (double *, int)
// my_rank is rank of process where this function is called
// p is # of processes
// com is a communicator for MPI

#include "mpi.h"
#include <stdio.h>

extern double MonteCarlo(double *a, double *b, int n, long int N, double (*fcn)(double *x, int n));

double parallelMonteCarlo(double *a, double *b, int n, long int N, double (*fcn)(double *x, int n), int my_rank, int p, MPI_Comm com)
{
	// integration summations
	double sum;
	double totalSum;

	// required variables for MPI
	int source;
	int dest = 0;
	int tag = 0;
	MPI_Status status;

	// If process != 0, perform monte carlo integration on N/p points
	// then send result to process 0
	if (my_rank > 0)
	{
		sum = MonteCarlo(a, b, n, (N/p), fcn);
		MPI_Send(&sum, 1, MPI_DOUBLE, dest, tag, com);
		// return sum;
	}

	// if process = 0, perform monte carlo integration on (N/P + N%p) points
	// then get results from other processes, add them as well
	else
	{
		totalSum = MonteCarlo(a, b, n, (N/p + N%p), fcn);
		//get results from other procs, add
		for(source = 1; source < p; source++)
		{
			MPI_Recv(&sum, 1, MPI_DOUBLE, source, tag, com, &status);
			totalSum += sum;
		}

		// Divide value by N to get avg result
		totalSum = totalSum/N;
		return totalSum;
	}
}