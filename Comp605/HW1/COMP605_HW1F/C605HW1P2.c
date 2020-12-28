/* Code for assignment 1 in Comp 605. 
This code contains a routine that utilizes the Jacobi iterative method to solve
the equation Ax = b where A is a matrix, x and b are vectors. The user can choose 
A to be a 3x3 or 4x4 matrix. The user also supplies the initial x values as a guess.

Author: Jon Parsons
Date: 2-10-19
*/

#include <stdio.h>
#include <stdlib.h>

int main(){
	
// Variable declerations
int Dime; 	// Dimensions of the matrices and vectors
int i;		// Looping integer
int j;		// Looping integer
int iternum;	// Determines the exit condition value, iteration is considered to have converged at this value of xvectn-xvect
int iter = 0; // Counts the iteration number
float err = 1.0; // Tracks the difference between xvectn and xvect	
float temp;	// Temporary values for iteration
	
//Begin user input
// Get Dimensions
printf("Input the dimensions to be used (3 or 4 recommended): \n");
scanf("%d", &Dime);
	
printf("Input the number of iterations desired \n");
scanf("%d", &iternum);

float **matA = malloc(Dime*sizeof(float*));
for (i = 0; i<Dime; i++) {
    matA[i] = malloc(Dime*sizeof(float));
}

float *xvect = malloc(Dime*sizeof(float));
float *xvectn = malloc(Dime*sizeof(float));
float *bvect = malloc(Dime*sizeof(float));
float *errvect = malloc(Dime*sizeof(float));

//Check for successful array allocation
if (!matA)
	{
	printf("Array A failed to allocate.\n");
	exit(0);
	}

if (!xvect)
	{
	printf("Array X failed to allocate.\n");
	exit(0);
	}

if (!xvect)
	{
	printf("Array Xn failed to allocate.\n");
	exit(0);
	}
	
if (!bvect)
	{
	printf("Array B failed to allocate.\n");
	exit(0);
	}

if (!errvect)
	{
	printf("Array Error failed to allocate.\n)");
	exit(0);
	}
	
FILE *FH = fopen("Jacerrs.dat","w");
	if (FH == NULL)
	{
		printf("File failed to open. \n");
		exit(0);
	}
	
fprintf(FH, "Iteration  Error \n");	
	
// Read in values of A
for (i=0; i<Dime; i++){
	for (j=0; j<Dime; j++){
		printf("Input value of matrix A at %d, %d\n",i,j);
		scanf("%f",&matA[i][j]);		
	}
}

// Read in values of x
for (i=0; i<Dime; i++){
		printf("Input initial value of X at %d\n",i);
		scanf("%f",&xvect[i]);
}
	
// Read in values of b
for (i=0; i<Dime; i++){
		printf("Input initial value of B at %d\n",i);
		scanf("%f",&bvect[i]);
}	
	
// Let user see the input
printf("\n Matrix A \n");
for (i=0; i<Dime; i++){
	for (j=0; j<Dime; j++){
		printf("%f  ", matA[i][j]);
	}
	printf("\n");
}

printf("\n");
	
printf("Initial values of X: \n");
for (i=0; i<Dime; i++){
		printf("%f  ", xvect[i]);
}

printf("\n");
	
printf("Initial values of B: \n");
for (i=0; i<Dime; i++){
		printf("%f  ", bvect[i]);
}

//Start of the while loop, exits when number of iterations reached.	
do{
	// Update iterations
	iter++; 
	printf("Iteration: %d \n ", iter);
	// Loops for the Jacobi method
	for (i=0; i<Dime; i++){
		temp = 0.0;
		for (j=0; j<Dime; j++){
			if (j != i){
				temp += matA[i][j]*xvect[j];
			}
		}
		xvectn[i] = (1.0/matA[i][i])*(bvect[i] - temp);
	}
	
	// Checks results
	for (i=0; i<Dime; i++){
		errvect[i] = 0.0;
		for (j=0; j<Dime; j++){
			errvect[i] += matA[i][j]*xvectn[j];
		}
	}
	
	// Output to file
	fprintf(FH, "%d, ( ", iter);
	for (i=0; i<Dime; i++){
		err += errvect[i];
		fprintf(FH, "%f, ",  xvectn[i]);
		xvect[i] = xvectn[i];
	}
	
			
}while(iter < iternum);

fprintf(FH, " )\n");
fprintf(FH," Results: (");
for (i=0; i<Dime; i++){
	fprintf(FH, "%f, ",  errvect[i]);
}
fprintf(FH, " ) \n");	
	
printf("Program Complete \n");
	
//Memory freeing
free(matA);
free(xvect);
free(xvectn);
free(errvect);
//End of Program
}

			
		
			
	
	
	
	
	
	

	
	
	
	
		
		
		
