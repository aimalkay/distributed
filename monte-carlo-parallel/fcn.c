// Aimal Khan SE4F03 Assignment 1 

// fcn.c

#include <math.h>
double fcn(double *x, int n)
{
	double pi = 3.14159265358979;
	double k = pi/2;
	double w = k;
	int i;
	for (i=0; i<n; i++)
		w *= x[i];
	double r = k*cos(w)-7*k*w*sin(w)-6*k*w*w*cos(w)+k*w*w*w*sin(w);
	return r;
}


