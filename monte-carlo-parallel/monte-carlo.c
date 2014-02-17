// Aimal Khan - 0946636 - SE4F03 Assignment 1 

// monte-carlo.c

// a is a pointer to n doubles, where a[i] stores a_i+1
// b is a pointer to n doubles, where a[i] stores b_i+1
// n is dimension n
// N is # of random points to be used in the integration
// fcn is a pointer to function returning double with arguments (double *, int)

#include <stdlib.h>
#include <stdio.h>
#define drand48() ((rand()/(double)(RAND_MAX)))

double MonteCarlo(double *a, double *b, int n, long int N, double (*fcn)(double *x, int n))
{
	double sum = 0;
	double *x = (double*)malloc(n*sizeof(double));
	
	// Iterate N times and add function value to sum
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < n; j++)
		{
			// random value between a and b, append to x
			x[j] = a[j] + drand48() * (b[j] - a[j]);
		}
		//call function for integration
		sum+=fcn(x, n);
	}

	free(x);
	return sum;
}