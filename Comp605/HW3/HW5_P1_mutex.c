/* Code for COMP 605 HW5, problem 1

Code will calculate pi from the integral of 
4/(1+x^2) on the bounds 0 - 1
in parallel using a pthread scheme. Mutex 
will be used as a resource control.

Author: Jon Parsons

compile using 
	gcc -o piThread_mutex.x HW5_P1_mutex.c -lm -lpthread
 */        

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

// Global variable declarations
const long NumDivs = 6400000; // Number of divisions in integration
const double dX = 1.0/6400000.0; // Step size along x 
long NumThreads;  // Number of threads, user input
double sum; // Global result of integration
pthread_mutex_t mutex; // pthread variabls

// Subroutine Declarations
void* SimpInt(void* rank);

double f(double x) // Called in SimpInt subroutine
/* Function to be evaluated. Returns pi when 
 integrated from 0 - 1    */
{	
	return 4.0/(1.0 + x*x);
}
/*---------------------------------------------------*/
int main() {
	long thread; // Local rank passed to integration subroutine
    pthread_t* thread_handles; // pthread variables
    double tDiff; // total time elapsed
	long tStart, tFin; // Initial and final times
	struct timeval timecheck; // Variable for timing routine
	
// User input
printf("Please input number of threads \n");
scanf("%ld", &NumThreads);

// pthread initialization 
thread_handles = (pthread_t*) malloc (NumThreads*sizeof(pthread_t)); 
   
// Initialize variables
pthread_mutex_init(&mutex, NULL);
sum = 0.0;

// Get start time in milliseconds
gettimeofday(&timecheck, NULL);
tStart = (long)timecheck.tv_sec*1000 + (long)timecheck.tv_usec/1000;
// Loop for calling simpint, calls same number of times as number of threads
for (thread = 0; thread < NumThreads; thread++)
{
	pthread_create(&thread_handles[thread], NULL,
	SimpInt, (void*)thread);  
}
// Return threads as well as final times
for (thread = 0; thread < NumThreads; thread++) 
{
	pthread_join(thread_handles[thread], NULL); 
	gettimeofday(&timecheck, NULL);
	tFin = (long)timecheck.tv_sec*1000 + (long)timecheck.tv_usec/1000;
	tDiff = tFin - tStart;
}
		
	// Output to screen 	
    printf("With Number of integral divisions = %ld terms,\n", NumDivs); 
    printf("pi found = %.15f, True = 3.14159265358979\n", sum);
    printf("Time Elapsed: %f milliseconds\n", tDiff);

	// Clean threads
	pthread_mutex_destroy(&mutex);
    free(thread_handles);
    return 0;
} 

/*------------------------------------------------------------------*/
void* SimpInt(void* rank) 
/* Function uses simpsons rule to find pi utilizing the 
	equation 4/(1+x^2) from the bounds 0 - 1
	Function is located in function f
*/
{
	long LocalRank = (long) rank; // cast rank to long
	long long i; // Looping integer
	long long LocalN = NumDivs/NumThreads; // Local number of divisions
	long long iStart = LocalN*LocalRank; // Local start of loop
	long long iFin = iStart + LocalN; // Local end of loop
	double LocalSum = 0.0; // Initialize the sum
	double xi; // Intermediate value of x

// Main loop
for (i = iStart; i < iFin; i++) 
{
	xi = i*dX; // Current value of x
	// Logic to determine weights for Simpson's rule
	if (i%2==0)
	{
		LocalSum += 2*f(xi); 
	}
	else
	{
		LocalSum += 4*f(xi);
	}
}
// Finalize local integration
LocalSum = (dX/3.0)*(f(0.0) + f(1.0) + LocalSum);
	// Exclude other threads until each is finished with their respective sections
    pthread_mutex_lock(&mutex);
	// After each thread is finished sum local results
    sum += LocalSum;
    // Remove locks
    pthread_mutex_unlock(&mutex);

   return NULL;
}

