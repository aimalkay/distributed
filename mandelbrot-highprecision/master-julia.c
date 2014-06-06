/*
 * -------------------------------------------------------------------------------------------------
 * Function: TaskMasterJulia
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
 * Outputs: int maxIterationCount - the maximum number of iterations required by any pixel in the
 *                                  process
 * -------------------------------------------------------------------------------------------------
 * This function designates process 0 as a master process. It is responsible for allocating work to
 * other processes and consolidating the work each slave process does when it returns to process 0.
 * Process 0 does not compute values for the image itself. When there are no more rows left to compute
 * in the image, the Master sends each slave process a DONE signal and exits.
 * 
 * All other processes are slaves. They are assigned work on a row-by-row basis. Slaves receive rows
 * as indexes and compute the Julia set for that row. It then passes back the row and waits for the 
 * next message from the Master. If there are more rows, the Master sends a new row index. If there
 * are no rows left, the Master sends a DONE message and the slave process exits.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpi.h>

#include "julia.h"

// Return types for Send/Receive
#define TYPEROW 0
#define TYPERETURN 1
#define TYPEDONE 2

// Define process 0 as master
#define MASTER 0

// Size of all messages
#define SIZE 1

// Booleans
#define FALSE 0
#define TRUE 1

long int TaskMasterJulia(mpf_t xmin, mpf_t xmax, unsigned long int xres, mpf_t ymin, mpf_t ymax, unsigned long int yres, mpf_t cr, mpf_t ci, 
	  int flag, int maxIterations, int *iterations, int my_rank, int p, MPI_Comm comm)
{
  int totalCount = 0;
  int done = FALSE;
  MPI_Status status;

  // Row number to pass to slave processes
  int *row;
  row = ( int* )malloc( sizeof(int) );
  assert(row != NULL);

  // Block for passing rows between processes
  int *block;
  block = ( int* )malloc( sizeof(int) * xres );
  assert(block != NULL);

  // Master process is only responsible for row allocation - does no work on Julia
  if (my_rank == MASTER)
  {   
    printf("Master process %d, ready to crack the whip! Allocating work...\n", my_rank);

    // FOR loop counter
    int i;
    
    // Used to calculate location in image (Big endian)
    int location;

    // Count how many rows each process completed
    int *processRows;
    processRows = ( int* )malloc( sizeof(int) * p );
    assert(processRows != NULL);

    // Keep track of which process has which row
    int *tracker;
    tracker = ( int* )malloc( sizeof(int) * p );
    assert(tracker != NULL);

    // Track image completion
    int sent = 0;
    int recv = 0;
    tracker[MASTER] = -1;
  
    // Initialize row counters for each process
    for (i = 0; i < p; i++) processRows[i] = 0;

    // Initially send one row to each slave process
    for (i = 1; i < p; i++)
    {
       *row = sent;
       MPI_Send(row, SIZE, MPI_INT, i, TYPEROW, comm);
       tracker[i] = *row;
       sent++;
    }
    
    // Have not heard about every row completion
    while (done == FALSE)
    {
       // Receive message from any process
       MPI_Recv(block, xres, MPI_INT, MPI_ANY_SOURCE, TYPERETURN, comm, &status);       

       // Make sure row is in bounds
       if (tracker[status.MPI_SOURCE] < yres)
       {
         // Update processed and received counters
         processRows[status.MPI_SOURCE]++;
         recv++;
         printf("\rCompleted: %lf%%", ((double)recv/yres)*100);

         // Put row data into image memory block
         location = tracker[status.MPI_SOURCE]*xres;
         for(i = 0; i < xres; i++) iterations[location + i] = block[i];

         // Received all rows from slave processes; send out DONE signal and exit
         if(recv == yres) 
         {
           done = TRUE;
           for(i = 1; i < p; i++) MPI_Send(row, SIZE, MPI_INT, i, TYPEDONE, comm);
         }
       
         // Have not sent all rows yet
         //if(done == FALSE)
         else
         {
           // Get next row, send to slave process, and update tracker and global sent counter
           *row = sent;
           MPI_Send(row, SIZE, MPI_INT, status.MPI_SOURCE, TYPEROW, comm);
           tracker[status.MPI_SOURCE] = *row;
	   sent++;
         }
       }
    }

    // Output how many rows each process completed
    for (i = 0; i < p; i++) printf("Rows completed on process %d: %d\n", i, processRows[i]);

    // Free memory on MASTER
    free(processRows);
    free(tracker);
  }
 
  // Slave process work on Julia one row at a time; Send message to master process when finished
  // Slaves are allocated additional work if there are remaining rows in Julia
  else
  {    
    printf("Slave process %d, reporting for duty!\n", my_rank);

    // Store iterations on a row
    int count;

    // Still rows to process
    while (done == FALSE)
    {
      // Get message from master process
      MPI_Recv(row, SIZE, MPI_INT, MASTER, MPI_ANY_TAG, comm, &status);
    
      if(status.MPI_TAG != TYPEDONE)
      {
        // Run Julia function, return block of iteration values
        count = julia(xmin, xmax, xres, xres, 0, ymin, ymax, SIZE, yres, *row, cr, ci, flag, maxIterations, block);
        totalCount += count;

        MPI_Send(block, xres, MPI_INT, MASTER, TYPERETURN, comm);
      }
      // Received DONE signal from MASTER - no more tasks
      else done = TRUE;
    }
  }

  // Free ALL OF THE MEMORY
  free(row);
  free(block);

  return totalCount;
}
