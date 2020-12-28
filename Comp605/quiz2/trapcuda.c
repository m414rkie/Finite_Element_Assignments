/* File:    trap.c
 * Purpose: Calculate definite integral using trapezoidal 
 *          rule.
 *
 * Input:   a, b, n
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Compile: gcc -g -Wall -o trap trap.c -lm
 * Usage:   ./trap
 *
 * Note:    The function f(x) is hardwired.
 *
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const int threadsPerBlock = 10;
const int blocksPerGrid = 10;

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x

*/
__device__ double f(double x) {
   double return_val;

   return_val = 4/(1+x*x);
   return return_val;
}  /* f */


double Trap(double a, double b, int n, double h);

int main(void) {
   double  integral;   /* Store result in integral   */
   double  a, b;       /* Left and right endpoints   */
   int     n;          /* Number of trapezoids       */
   double  h;          /* Height of trapezoids       */
   double  err;
   
   a = 0;
   b = 1;
   n = 100000;

   h = (b-a)/n;
   integral = Trap<<<blocksPerGrid,threadsPerBlock>>>(a, b, n, h);
   err = fabs(integral - M_PI);
   
   printf("With n = %d trapezoids, our estimate\n", n);
   printf("of the integral is %.15f\n",integral);
   
   printf("Numerical error is %.15f\n", err);

   return 0;
}  /* main */

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids
 * Input args:  a, b, n, h
 * Return val:  Estimate of the integral 
 */
double Trap(double a, double b, int n, double h) {
   double integral;
    __shared__ float cache[threadsPerBlock];
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;

   integral = (f(a) + f(b))/2.0;
	
	while(tid < n)
	{
		integral += f(a+tid*h);
	}

   integral = integral*h;	
	
    // set the cache values
    cache[cacheIndex] = integral;
    
    // synchronize threads in this block
    __syncthreads();

    // for reductions, threadsPerBlock must be a power of 2
    // because of the following code
    int i = blockDim.x/2;
    while (i != 0) {
        if (cacheIndex < i)
            cache[cacheIndex] += cache[cacheIndex + i];
        __syncthreads();
        i /= 2;
    }

   return integral;
}  /* Trap */

