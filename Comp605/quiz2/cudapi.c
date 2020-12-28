/* Code for COMP 605 HW5, problem 1

Code will calculate pi from the integral of 
4/(1+x^2) on the bounds 0 - 1
using CUDA methodology.

Author: Jon Parsons

compile using 
	nvcc -o CudaPI.x  cudapi.c
*/        

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int threadsPerBlock = 10;
const int blocksPerGrid = 10;

__device__ double f(double x) // Called in SimpInt subroutine
/* Function to be evaluated. Returns pi when 
 integrated from 0 - 1    */
{	
	return 4.0/(1.0 + x*x);
}

__global__ __device__ void int_kernel(double *f_vals, double b, double h, int n)
{
	int tidx = blockIdx.x * blockDim.x + threadIdx.x;
	double x_loc = b + (float) tidx*h;
		
	if (idx<n)
	{
		f_vals[tidx] = f(x) + f(x+h);
	}
	
}

__host__ double TrapInt(double a, double b, int n) 
/* Function uses simpsons rule to find pi utilizing the 
	equation 4/(1+x^2) from the bounds 0 - 1
	Function is located in function f
*/
{
	double h = (b-a)/n;
	
	double* f_vals = (float *) malloc(n*sizeof(double));
	
	cudaMalloc((void **)&f_vals_loc,(n*sizeof(double)));
	
	int_kernel<<<threadsPerBlock,blocksPerGrid>>>(f_vals_loc, a, h, n);
	
	cudaMemcpy(f_vals,f_vals_loc, sizeof(float)*n, cudaMemcpyDeviceToHost);
	
	double sum=0.0:
	for (int i=0; i<n; i++) sum += f_vals[i];
	
	sum *= h*0.5;
	
	free(f_vals);
	cudaFree(f_vals_loc);
	
	return sum;
	
}

/*---------------------------------------------------*/
int main() {

	double integral = TrapInt(a, b, n);
	
	printf ("Integral result %f", integral);
	
    return 0;
} 

/*------------------------------------------------------------------*/


