/*
 * -------------------------------------------------------------------------------------------------
 * Function: parallelJulia
 * Inputs: mpf_t xmin, xmax - x coordinates
 *         unsigned long int xres - the width of the complete image
 *         mpf_t ymin, ymax - y coordinates
 *	   int yblock - the height of the block julia is computing values for 
 *         unsigned long int yres - the height of the complete image
 *         int starty - y offset of the memory block that julia is working on
 *         mpf_t cr, ci - values of the imaginary number c + ci
 *         int flag - indicates if the image is Mandelbrot or Julia set
 *         int maxIterations - maximum number of hops to try to exit the unit circle
 *         int *iterations - the memory block that julia is working on
 *         int my_rank - the id of the current process
 *         int p - the total number of processes running
 *         MPI_Comm comm - the MPI communicator of the program
 * Outputs: int maxCount - the maximum number of iterations required by any pixel in the
 *                         process
 * -------------------------------------------------------------------------------------------------
 * This function determines which Julia computation algorithm to use based on the number of 
 * processes:
 *  - # Processes = 1: Serial program; call julia function directly
 *  - # Processes = 2: Not enough processes to require a task master; send to BlockPartitionJulia
 *  - # Processes > 2: Enough processes to require a task master; send to TaskMasterJulia
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpi.h>

#include "julia.h"

long int parallelJulia(mpf_t xmin, mpf_t xmax, unsigned long int xres, mpf_t ymin, mpf_t ymax, unsigned long int yres, mpf_t cr, mpf_t ci, 
	  int flag, int maxIterations, int *iterations, int my_rank, int p, MPI_Comm comm)
{
  long int count;

  if (p == 1)
  {
    if(my_rank == 0) printf("Single process - serial version\n\n");
    printf("Process %d...aren't you happy I'm here?\n", my_rank);

    count = julia(xmin, xmax, xres, xres, 0, ymin, ymax, yres, yres, 0, cr, ci, flag, maxIterations, iterations);
  }
  else if (p == 2)
  {
    if(my_rank == 0) printf("Not enough processes - divide image at middle height and use scatterv/gatherv\n\n");
    count = BlockPartitionJulia(xmin, xmax, xres, ymin, ymax, yres, cr, ci, flag, maxIterations, iterations, my_rank, p, comm);
  }
  else
  {
    if(my_rank == 0) printf("Sufficient processes - run process 0 as task master\n\n");
    count = TaskMasterJulia(xmin, xmax, xres, ymin, ymax, yres, cr, ci, flag, maxIterations, iterations, my_rank, p, comm);
  }

  return count;
}
