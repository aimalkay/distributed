/*
 * -------------------------------------------------------------------------------------------------
 * Function: getParams
 * Inputs: char **argv - a list of arguements from the command line
 *	   int *flag - a pointer to an INT in memory to store a set flag
 *         mpf_t *cr, *ci - a pointer to a pair of values representing a complex number cr + ci
 *         mpf_t *x, *xr - a pointer to the center and radius of the x circle
 *         mpf_t *y, *yr - a pointer to the center and radius of the y circle
 *         unsigned long int *height - a pointer to an INT to store the width of the image
 *         unsigned long int *width - a pointer to an INT to store the height of the image
 *         int *maxiter - a pointer to an INT to store the maximum number of iterations to use
 *         char **image - a pointer to a file name to store the image pixels; .bmp extension
 * -------------------------------------------------------------------------------------------------
 * This function takes in an array storing arguements from a file and parses them into the
 * program's memory. There are no values to return because all values are returned by reference.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#define SIZE 100

void getParams(char **argv, int *flag, mpf_t *cr, mpf_t *ci, mpf_t *x, mpf_t *y, mpf_t *xr, mpf_t *yr, unsigned long int *height, unsigned long int *width, int *maxiter, char **image)
{
  char data[SIZE];
  char filename[SIZE];

  FILE *params;
  params = fopen(argv[1], "r");

  if (params != NULL)
  {
    if (fgets(data, SIZE, params) != NULL) *flag = strtol(data, NULL, 0);
    if (fgets(data, SIZE, params) != NULL) mpf_set_str(*cr, data, 10);
    if (fgets(data, SIZE, params) != NULL) mpf_set_str(*ci, data, 10);
    if (fgets(data, SIZE, params) != NULL) mpf_set_str(*x, data, 10);
    if (fgets(data, SIZE, params) != NULL) mpf_set_str(*y, data, 10);
    if (fgets(data, SIZE, params) != NULL) mpf_set_str(*xr, data, 10);
    if (fgets(data, SIZE, params) != NULL) mpf_set_str(*yr, data, 10);
    if (fgets(data, SIZE, params) != NULL) *height = strtol(data, NULL, 0);
    if (fgets(data, SIZE, params) != NULL) *width = strtol(data, NULL, 0);
    if (fgets(data, SIZE, params) != NULL) *maxiter = strtol(data, NULL, 0);
    if (fgets(data, SIZE, params) != NULL) 
    {
      sscanf(data,"%s", filename);
      *image = malloc(strlen(filename));
      strcpy(*image, filename);
    }

    fclose(params);
  }
  else perror("Error opening file\n");

  return;
}
