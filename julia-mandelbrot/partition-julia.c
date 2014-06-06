/*
 * -------------------------------------------------------------------------------------------------
 * Function: BlockPartitionJulia
 * Inputs: mpf_t xmin, xmax - x coordinates
 *         unsigned long int xres - the width of the complete image
 *         mpf_t ymin, ymax - y coordinates
 *	   unsigned long int yres - the height of the complete image
 *         mpf_t cr, ci - values of the imaginary number cr + ci
 *         int flag - indicates if the image is Mandelbrot or Julia set
 *         int maxIterations - maximum number of hops to try to exit the unit circle
 *         int *iterations - the memory block that julia is working on
 *         int my_rank - the id of the current process
 *         int p - the total number of processes running
 *         MPI_Comm comm - the MPI communicator of the program
 * Outputs: int maxCount - the maximum number of iterations required by any pixel in the
 *                         process
 * -------------------------------------------------------------------------------------------------
 * This function evenly divides the image to calculate into blocks and sends them to each process
 * using Scatterv. The image blocks are recombined on process 0 using Gatherv.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpi.h>

#include "julia.h"

long int BlockPartitionJulia(mpf_t xmin, mpf_t xmax, unsigned long int xres, mpf_t ymin, mpf_t ymax, unsigned long int yres, mpf_t cr, mpf_t ci, int flag, int maxIterations, int *iterations, int my_rank, int p, MPI_Comm comm)
{
  int i;

  int *block_size;
  block_size = ( int* )malloc( sizeof(int) * p );
  assert(block_size != NULL);

  int *sendElements;
  sendElements = ( int* )malloc( sizeof(int) * p );
  assert(sendElements != NULL);

  int remaining = yres % p;

  // For Scatterv and Gatherv
  int *displacement;
  displacement = (int*)malloc( sizeof(int) * p );
  assert(displacement != NULL);

  // For julia.c
  int *offset;
  offset = ( int* )malloc( sizeof(int) * p );
  assert(offset != NULL);

  printf("Process %d reporting for duty!\n", my_rank);

  for (i = 0; i < p; i++)
  {
    // Determine block size
    block_size[i] = yres / p;
    if (i < remaining) block_size[i]++;

    // Determine how many elements to send to each process
    sendElements[i] = block_size[i] * xres;

    // Determine displacement from iterations[0, 0] for Scatterv and Gatherv
    if (i == 0) displacement[i] = 0;
    else displacement[i] = displacement[i - 1] + sendElements[i - 1];

    // Determine row offset for julia.c
    if (i == 0) offset[i] = 0;
    else offset[i] = offset[i - 1] + block_size[i];
  }

  if(my_rank == 0)
    for (i = 0; i < p; i++) 
    {
	printf("Process %d: Elements = %d, Displacement = %d \n", i, sendElements[i], displacement[i]);
	printf("Process %d: Block Size = %d, Offset = %d \n", i, block_size[i], offset[i]);
    }

  // Allocate space for local arrays
  int *block;
  block = ( int* )malloc( sizeof(int) * sendElements[my_rank] );
  assert(block != NULL);

  // Send data to processes
  MPI_Scatterv(iterations, sendElements, displacement, MPI_INT, block, sendElements[my_rank], MPI_INT, 0, comm);  

  // Run julia
  int xblock = xres;
  long int count = julia(xmin, xmax, xblock, xres, 0, ymin, ymax, block_size[my_rank], yres, offset[my_rank], cr, ci, flag, maxIterations, block);

  // Gather blocks back into interations
  MPI_Gatherv(block, sendElements[my_rank], MPI_INT, iterations, sendElements, displacement, MPI_INT, 0, comm);

  // Free ALL OF THE MEMORY!!!
  free(block_size);
  free(sendElements);
  free(displacement);
  free(offset);
  free(block);

  return count;
}
