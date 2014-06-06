/*
 * -------------------------------------------------------------------------------------------------
 * This header file contains function signatures for each function that is required to successfully
 * run the Julia program.
*/

long int TaskMasterJulia(mpf_t xmin, mpf_t xmax, unsigned long int xres, mpf_t ymin, mpf_t ymax, unsigned long int yres, mpf_t cr, mpf_t ci, 
	  int flag, int maxIterations, int *iterations, int my_rank, int p, MPI_Comm comm);

long int BlockPartitionJulia(mpf_t xmin, mpf_t xmax, unsigned long int xres, mpf_t ymin, mpf_t ymax, unsigned long int yres, mpf_t cr, mpf_t ci, int flag, int maxIterations, int *iterations, int my_rank, int p, MPI_Comm comm);

long int parallelJulia(mpf_t xmin, mpf_t xmax, unsigned long int xres, mpf_t ymin, mpf_t ymax, unsigned long int yres, mpf_t cr, mpf_t ci, 
	  int flag, int maxIterations, int *iterations, int my_rank, int p, MPI_Comm comm);

long int julia(mpf_t xmin, mpf_t xmax, int xblock, unsigned long int xres, int startx, mpf_t ymin, mpf_t ymax, int yblock, unsigned long int yres, int starty, mpf_t cr, mpf_t ci, int flag, int maxIterations, int *iterations);

void getParams(char **argv, int *flag, mpf_t *cr, mpf_t *ci, mpf_t *x, mpf_t *y, mpf_t *xr, mpf_t *yr, unsigned long int *height, unsigned long int *width, int *maxiter, char **image);

void saveBMP(char* filename, int* result, int width, int height);
