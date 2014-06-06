/*
 * -------------------------------------------------------------------------------------------------
 * Function: main
 * Inputs: int argc - the number of arguements passed in from the command line
 *         char *argv - the list of arguements passed in from the command line
 * -------------------------------------------------------------------------------------------------
 * This function initialize memory blocks that will be used for the duration of the program. It then
 * calls getParams to parse the command line arguements into the allocated memory before initializing
 * the MPI environment.
 * 
 * After MPI is initialize, a timer is started before the processes begin their Julia set 
 * calculations. When each process finishes, the timer is stopped and the statistics are collected on
 * process 0 for output to a stats file. Process 0 is also responsible for converting the iterations
 * calculated by Julia and converting them into .bmp files.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpi.h>
#include <string.h>

#include "julia.h"

int main(int argc, char *argv[])
{
  int maxiter, flag;
  unsigned long int width, height;
  char *image;
  //long int precision, temp;

  mpf_t cr, ci, x, y, xr, yr, xmin, xmax, ymin, ymax;
  mpf_set_default_prec(300);
  //mpf_inits(cr, ci, x, y, xr, yr, xmin, xmax, ymin, ymax, (mpf_t *) 0);
  mpf_init(cr);
  mpf_init(ci);
  mpf_init(x);
  mpf_init(y);
  mpf_init(xr);
  mpf_init(yr);
  mpf_init(xmin);
  mpf_init(xmax);
  mpf_init(ymin);
  mpf_init(ymax);

  int comm_sz, my_rank;
  double t1, t2, delta, maxTime;
  long int totalIterations;

  // Get and parse the program parameters
  getParams(argv, &flag, &cr, &ci, &x, &y, &xr, &yr, &width, &height, &maxiter, &image);

  // xmin and xmax
  mpf_sub(xmin, x, xr);
  mpf_add(xmax, x, xr);

  // ymin and ymax
  mpf_sub(ymin, y, yr);
  mpf_add(ymax, y, yr);
  
  // Allocate space for the image
  int *iterations = (int*)malloc( sizeof(int) * width * height );
  assert(iterations != NULL);

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0)
  {
    int n = 30;

    printf("Flag = %d\n", flag);

    gmp_printf ("cr = %+.*Ff\t", n, cr);
    gmp_printf ("ci = %+.*Ff\n", n, ci);

    gmp_printf ("x = %+.*Ff\t", n, x);
    gmp_printf ("xr = %+.*Ff\n", n, xr);
    gmp_printf ("xmin = % .*Ff\t", n, xmin);
    gmp_printf ("xmax = % .*Ff\n", n, xmax);

    gmp_printf ("y = %+.*Ff\t", n, y);
    gmp_printf ("yr = %+.*Ff\n", n, yr);
    gmp_printf ("ymin = % .*Ff\t", n, ymin);
    gmp_printf ("ymax = % .*Ff\n", n, ymax);

    printf("Height = %ld\t", height);
    printf("Width = %ld\t", width);
    printf("maxiter = %d\t", maxiter);
    printf("Image = %s\n", image);
  }

  printf("Process %d waiting for other processes\n", my_rank);
  MPI_Barrier(MPI_COMM_WORLD);

  t1 = MPI_Wtime();

  /* Compute Julia set */
  long int count;
  count = parallelJulia(xmin, xmax, width, ymin, ymax, height, cr, ci, flag, maxiter, iterations, my_rank, comm_sz, MPI_COMM_WORLD);

  t2 = MPI_Wtime();

  printf("Process %d waiting for Julia set completion\n", my_rank);
  MPI_Barrier(MPI_COMM_WORLD);
  delta = t2 - t1;

  MPI_Reduce(&count, &totalIterations, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&delta, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (my_rank == 0)
  {
    /* save our picture for the viewer */
    printf("\nMaster process %d creating image...\n", my_rank);
    saveBMP(image, iterations, width, height);
    printf("%d  %lf  %ld\n", comm_sz, maxTime, totalIterations);
  }

  MPI_Finalize();

  // Free reserved memory
  //mpf_clears(cr, ci, x, y, xr, yr, xmin, xmax, ymin, ymax, (mpf_t *) 0);
  mpf_clear(cr);
  mpf_clear(ci);
  mpf_clear(x);
  mpf_clear(y);
  mpf_clear(xr);
  mpf_clear(yr);
  mpf_clear(xmin);
  mpf_clear(xmax);
  mpf_clear(ymin);
  mpf_clear(ymax);
  free(iterations);

  return 0;
}

 /*
 
  //Set xmin and xmax
  precision = mpf_get_prec(x);
  temp = mpf_get_prec(xr);
  if (temp > precision)
  {
	precision = temp;
  }
  
  mpf_init2(xmin, precision);
  mpf_init2(xmax, precision);
  
  mpf_add(xmax, x, xr);
  mpf_sub(xmin, x, xr);
  
  //Set ymin and ymax
  precision = mpf_get_prec(y);
  temp = mpf_get_prec(yr);
  if (temp > precision)
  {
	precision = temp;
  }
  
  mpf_init2(ymin, precision);
  mpf_init2(ymax, precision);
  
  mpf_set_default_prec(precision);
  
  mpf_add(ymax, y, yr);
  mpf_sub(ymin, y, yr);
  */
