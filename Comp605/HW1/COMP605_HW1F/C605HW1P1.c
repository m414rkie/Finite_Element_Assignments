/* Code for assignment 1 in Comp 605. 
This code contains a routine that multiplies two matrices together and outputs to 
console. Arrays are dynamically sized from user input. Values are user input as well.

Author: Jon Parsons
Date: 2-10-19
*/

#include <stdio.h>
#include <stdlib.h>

int main(){
	
// Variable declerations
int DimA; 	// Dimension 1 of Matrix A
int Dim2;	// Second Dimension of Matrix A and first dimension of matrix B
int DimB;	// Second Dimension of Matrix B
int i;		// Looping integer
int j;		// Looping integer
int k;		// Looping integer

//Begin user input
// Get Dimensions
printf("Input the 1st dimension of matrix A: \n");
scanf("%d", &DimA);
printf("Input the 2nd dimension of matrix A ( Will be 1st dimension of matrix B) \n");
scanf("%d", &Dim2);
printf("Input the 2nd Dimension of matrix B: \n");
scanf("%d", &DimB);

float **matA = malloc(DimA*sizeof(float*));
for (i = 0; i<Dim2; i++) {
    matA[i] = malloc(Dim2*sizeof(float));
}

float **matB = malloc(Dim2*sizeof(float*));
for (i = 0; i<Dim2; i++) {
    matB[i] = malloc(DimB*sizeof(float));
}

float **matC = malloc(DimA*sizeof(float*));
for (i = 0; i<DimB; i++) {
    matC[i] = malloc(DimB*sizeof(float));
}

//Check for successful array allocation
if (!matA)
	{
	printf("Array A failed to allocate.\n");
	exit(0);
	}

if (!matB)
	{
	printf("Array B failed to allocate.\n");
	exit(0);
	}

if (!matC)
	{
	printf("Array C failed to allocate.\n");
	exit(0);
	}

// Initialize C
for (i=0; i<DimA; i++){
	for (j=0; j<DimB; j++){
		matC[i][j] = 0.0;
	}
}

// Read in values of A
for (i=0; i < DimA; i++){
	for (j=0; j <Dim2; j++){
		printf("Input value of matrix A at %d, %d\n",i,j);
		scanf("%f",&matA[i][j]);		
	}
}

// Read in values of B
for (i=0; i < Dim2; i++){
	for (j=0; j <DimB; j++){
		printf("Input value of matrix B at %d, %d\n",i,j);
		scanf("%f",&matB[i][j]);		
	}
}

// Let user see the input
printf("\n Matrix A\n");
for (i=0; i < DimA; i++){
	for (j=0; j <Dim2; j++){
		printf("%f  ", matA[i][j]);
	}
	printf("\n");
}

printf("\n");
	
printf("Matrix B\n");
for (i=0; i < Dim2; i++){
	for (j=0; j <DimB; j++){
		printf("%f  ", matB[i][j]);
	}
	printf("\n");
}
	
// Begin Matrix Multiplication Part
// Simple matrix multiplication algorithm

for (i=0; i<DimA; i++){
	for (j=0; j<DimB; j++){
		for (k=0; k<Dim2; k++){
			matC[i][j] += matA[i][k]*matB[k][j];
		}
	}
}
	
printf("\n Calculation Complete \n");
// Output
printf("Results:\n");
for (i=0; i < DimA; i++){
	for (j=0; j < DimB; j++){
		printf("%f  ", matC[i][j]);
	}
	printf("\n");
}

printf("Program Complete \n");
	
//Memory freeing
free(matA);
free(matB);
free(matC);

//End of Program
}

			
		
			
	
	
	
	
	
	

	
	
	
	
		
		
		
