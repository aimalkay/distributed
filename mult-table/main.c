// SE4F03 Final Project
// Aimal Khan
// Sean McLellan

// main.c

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include "mpi.h"



int main (int argc, char *argv[])
{
	//MPI Initialization
	int p, my_rank;
	MPI_Status status;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//Dimension of table
	long long int n;
	//Processor specific unique count
	long long int uniqueCount = 0;
	//Total unique count
	long long int totalUniqueCount = 0;
	int row,col,i;
	//Determines whether final row of segment has been reached
	bool bEndOfSegment = false;
	
	/* SET THESE BEFORE COMPILE */
	int ranks_per_node = p; //if only 1 node, just number of cores

	//Memory use in MB
	long long int MEMORY = 2 * 1024; //Memory per node
	long long int BUFFER = 20; //Amount of memory to afford buffer per rank (in MB)

	//Retrieve parameters
	if (my_rank == 0)
	{
		if (argc == 2)
		{
			if (argv[1] != NULL)
			{
				n = atol(argv[1]);
			}
		}
	}

	//Broadcast n to all processors
	MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	//Only adjust if you know how much memory your cores will have available to them.
	//If this is set to a value, MEMORY and BUFFER are ignored.
	long long int MAX_SIZE = 26000000;

	//(((MEMORY * 1024 * 1024)/ranks_per_node) - (n * 2 * sizeof(long long int)) - BUFFER * 1024 * 1024)/(sizeof(long long int));
	
	//printf("MAX ARRAY SIZE: %lld\n", MAX_SIZE);

	//The limit on memory must be at least as large as n
	//Increasing MAX_SIZE increases memory usage, but increases performance
	//This adjustment may cause unexpected behaviour. Set the memory before compile.
	if (MAX_SIZE < n)
		MAX_SIZE = n;
	
	long long int *segments;
	segments = (long long int*)malloc(2 * n * sizeof(long long int));
	memset (segments, 0, 2 * n * sizeof(long long int));
	
	//Distribute workload evenly across processors
	//Each processor calculates own segment
	int processorsUsed = p;
	if (n < p)
		processorsUsed = n;

	//Set MAX_SIZE
	if (MAX_SIZE > (n * n))
	{
		MAX_SIZE = (n * n) / processorsUsed + (n * n) % processorsUsed;
	}
	
	//The row that the last segment ends on, this checks termination
	long long int prevRowEnd = 0;
	//Total elements given to segment
	long long int totalElements = 0;
	int processor;
	
	long long int r = 0;
	long long int iteration = 0;
	long long int highestIteration = 0;
	
	while (prevRowEnd < n)
	{
		//printf("%lld\n", prevRowEnd);
		for (processor = 0; processor < processorsUsed; processor++)
		{
			//Reset each iteration
			totalElements = 0;
			//Set initial start row for segment, this is the min
			if (my_rank == processor)
			{
				segments[iteration] = prevRowEnd;
				//printf("rank: %d, start:%lld\n", my_rank, segments[iteration]);
			}
			
			for (r = prevRowEnd; r < n; r++)
			{
				//The step is the size of the n value
				if (totalElements + n <= MAX_SIZE)
				{
					totalElements += n;
				}
				else
				{
					//If the total number of elements stored in segment is too large to add another row
					//break the loop, set the prevEndRow to r.
					break;
				}
			}

			prevRowEnd = r;
			if (my_rank == processor)
			{
				//This is the max
				segments[iteration + 1] = prevRowEnd;
				//printf("rank: %d, end:%lld\n", my_rank, segments[iteration + 1]);
				//Track the highest iteration for this processor
				highestIteration = iteration + 1;
			}
		}
		//Increase the iteration by 2, each iteration stores a min and max
		iteration+=2;
	}

	//CPU Time
	double time;
	clock_t start, end;
	start = clock();

	//Iterate through all segments previously defined above unique to each processor.
	for (i = 0; i <= highestIteration; i+=2)
	{
		//Tracking the appearance of all products in the segment
		long long int *productCount;

		//Reset after each segment
		bEndOfSegment = false;

		long long int startRow, endRow;
		startRow = segments[i];
		endRow = segments[i+1];
		
		long long int max;
		long long int min;

		max = endRow * n;
		min = startRow * n;

		long long int range;
		range = max - min;

		//Buffer for segment is at most of size n for each segment
		productCount = (long long int*)malloc(range * sizeof(long long int));
		memset (productCount, 0, range * sizeof(long long int));

		for (row = startRow + 1; row <= n; row++)
		{
			for (col = row; col <= n; col++)
			{
				long long int product = row * col;
				//printf("ROW:%lld COL:%lld P:%lld MAX:%lld\n", row, col, product, max);

				if (product > min && product <= max)
				{
					//Gives each value in segment a unique value based on max value
					//of segment for tracking count
					if (productCount[max - product] <= 0)
					{
						//printf("%lld, %lld, %lld\n", product, min, max);
						uniqueCount++;
						productCount[max - product]++;
					}
				}
		
				if (product > max)
				{
					//printf("%lld, %lld, %lld\n", row, col, max);
					if (row == col)
						bEndOfSegment = true;
					break;
				}
			}

			if (bEndOfSegment)
				break;
		}

		//Release memory
		free (productCount);
	}

	//Stop clock 
	end = clock();
	time = ((double)(end - start)) / CLOCKS_PER_SEC;
	
	double maxTimePerCPU = 0;

	MPI_Reduce(&uniqueCount, &totalUniqueCount, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &maxTimePerCPU, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (my_rank == 0)
	{
		printf("Time: %f\n", maxTimePerCPU);
		printf("M(N): %lld\n", totalUniqueCount);
		printf("M(N)/N^2: %f\n", (float)totalUniqueCount/(float)(n*n));
	}

	MPI_Finalize();
	return 0;
}
