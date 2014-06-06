/*
 * -------------------------------------------------------------------------------------------------
 * Function: julia
 * Inputs: mpf_t xmin, xmax - x coordinates
 *	   int xblock - the width of the block julia is computing values for 
 *         unsigned long int xres - the width of the complete image
 *         int startx - x offset of the memory block that julia is working on
 *         mpf_t ymin, ymax - y coordinates
 *	   int yblock - the height of the block julia is computing values for 
 *         unsigned long int yres - the height of the complete image
 *         int starty - y offset of the memory block that julia is working on
 *         mpf_t cr, ci - values of the imaginary number c + ci
 *         int flag - indicates if the image is Mandelbrot or Julia set
 *         int maxIterations - maximum number of hops to try to exit the unit circle
 *         int *iterations - the memory block that julia is working on
 * Outputs: int maxIterationCount - the maximum number of iterations required by any pixel in the
 *                                  memory block
 * -------------------------------------------------------------------------------------------------
 * This function generates iteration values for pixels in an image of resolution xres x yres. The
 * memory block iterations is a xblock x yblock section of the complete image. The values startx and
 * starty tell julia where the block is located in the image.
 * Julia tried maxIterations times to leave the unit circle for each pixel and records that value in
 * iterations. The memory block iterations does not need to be explicitly returned because it is 
 * passed by reference.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpi.h>

#include "julia.h"

long int julia(mpf_t xmin, mpf_t xmax, int xblock, unsigned long int xres, int startx, mpf_t ymin, mpf_t ymax, int yblock, unsigned long int yres, int starty, mpf_t cr, mpf_t ci, int flag, int maxIterations, int *iterations)
{
  /* Maximum radius of the unit circle */
  const double maxRadius = 2.0;
  
  /* Set precision of numbers */
  long int precision;
  int compare;
  precision = mpf_get_prec(xmax);
  mpf_set_default_prec(precision);

  /* Counters */
  int iteration;
  long int iterationCount = 0;
  int i, j;
  
  /* Distance variables */
  mpf_t xgap, ygap;

  //mpf_inits(xgap, ygap, (mpf_t *) 0);
  mpf_init(xgap);
  mpf_init(ygap);

  /* Complex calculation variables */
  mpf_t zinitReal, zinitImag;
  mpf_t z0Real, z0Imag;
  mpf_t zReal, zImag;
  mpf_t cReal, cImag;
  mpf_t tempReal, tempImag;
  mpf_t magnitude;
   
  //mpf_inits(zinitReal, zinitImag, z0Real, z0Imag, zReal, zImag, cReal, cImag, tempReal, tempImag, magnitude, (mpf_t *) 0);
  mpf_init(zinitReal);
  mpf_init(zinitImag);
  mpf_init(z0Real);
  mpf_init(z0Imag);
  mpf_init(zReal);
  mpf_init(zImag);
  mpf_init(cReal);
  mpf_init(cImag);
  mpf_init(tempReal);
  mpf_init(tempImag);
  mpf_init(magnitude);
  
  /* Converting coordinate to complex space */
  mpf_sub(xgap, xmax, xmin);    // xgap = (x[1] - x[0]) / xres;
  mpf_div_ui(xgap, xgap, xres);
  
  mpf_sub(ygap, ymax, ymin);    // ygap = (y[1] - y[0]) / yres;
  mpf_div_ui(ygap, ygap, yres);
  
  /* Calculate Julia */
  for (j = 0; j < yblock; j++)
    { 
      for (i = 0; i < xblock; i++)    
	{
	  /* pixel to coordinates, base values for sets */
	  mpf_mul_ui(tempReal, xgap, (i+startx));  // zinitReal = x[0] + (i + startx) * xgap;
	  mpf_add(zinitReal, xmin, tempReal);
	  
	  mpf_mul_ui(tempImag, ygap, (j+starty));  // zinitImag = y[0] + (j + starty) * ygap;
	  mpf_add(zinitImag, ymin, tempImag);
    
          mpf_set(cReal, cr);                  // double complex C = c[0] + c[1]*I;
          mpf_set(cImag, ci);

	  /* if flag=0, z = z0, flag=1, z = C */
	  mpf_mul_ui(tempReal, cReal, flag);      //z0Real = flag*cReal + (1 - flag)*zinitReal;
	  mpf_mul_ui(tempImag, zinitReal, (1-flag));
	  mpf_add(z0Real, tempReal, tempImag);
		  
	  mpf_mul_ui(tempReal, cImag, flag);     //z0Imag = flag*cImag + (1 - flag)*zinitImag;
	  mpf_mul_ui(tempImag, zinitImag, (1-flag));
	  mpf_add(z0Imag, tempReal, tempImag);
	  
	  /* Determine how long it takes to leave the unit circle */
	  iteration = 0;
	  
          mpf_set(zReal, zinitReal);  // double complex z  = zinit;
          mpf_set(zImag, zinitImag);

          /* Calculate magnitude */
          mpf_mul(tempReal, zReal, zReal); // magnitude = (zReal*zReal) + (zImag*zImag);
	  mpf_mul(tempImag, zImag, zImag);
	  mpf_add(magnitude, tempReal, tempImag);
		  
          /* Check if the point is outside the unit circle */
	  compare = mpf_cmp_d(magnitude, (maxRadius*maxRadius));

	  while (compare < 0 && iteration < maxIterations)
	  {
            iteration++;

	    /* z = z*z+z0 = (a + bi)^2 = ([a^2 - b^2] + z0[real]) + ([2ab] + z0[imag])i */
              
            /* Real part - a^2 - b^2 */
	    mpf_mul(tempReal, zReal, zReal);       // tempReal = (zReal*zReal) - (zImag*zImag);
            mpf_mul(tempImag, zImag, zImag);
	    mpf_sub(tempReal, tempReal, tempImag);
			  
	    /* Imaginary part - 2ab */
            mpf_mul_ui(tempImag, zReal, 2);        // tempImag = 2 * zReal * zImag; 
	    mpf_mul(tempImag, tempImag, zImag);

            mpf_add(zReal, tempReal, z0Real);      // zReal = tempReal + z0Real;
            mpf_add(zImag, tempImag, z0Imag);      // zImag = tempImag + z0Imag;
			  
	    /* Calculate magnitude and check if the point is outside the unit circle */
	    mpf_mul(tempReal, zReal, zReal);          //magnitude = (zReal*zReal) + (zImag*zImag);
	    mpf_mul(tempImag, zImag, zImag);
	    mpf_add(magnitude, tempReal, tempImag);
	    compare = mpf_cmp_d(magnitude, (maxRadius*maxRadius));
	   }
	  
	  /* Count how many iterations are performed */
	  iterationCount += iteration;
	  
	  /* Calculate storage location and record iteration count for pixel */
	  int *p = iterations + j*xres+i;
	  *p = iteration;
	}
    }

  //mpf_clears(zinitReal, zinitImag, z0Real, z0Imag, zReal, zImag, cReal, cImag, tempReal, tempImag, magnitude, (mpf_t *) 0);
  mpf_clear(zinitReal);
  mpf_clear(zinitImag);
  mpf_clear(z0Real);
  mpf_clear(z0Imag);
  mpf_clear(zReal);
  mpf_clear(zImag);
  mpf_clear(cReal);
  mpf_clear(cImag);
  mpf_clear(tempReal);
  mpf_clear(tempImag);
  mpf_clear(magnitude);
  mpf_clear(xgap);
  mpf_clear(ygap);

  return iterationCount;
}
