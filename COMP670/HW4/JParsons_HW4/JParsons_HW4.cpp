/**
Program written for homework 4, COMP670.

This program will solve a 2D advection-diffusion system using the Mimetic
Operators Library Enhanced (MOLE)

Version 1.0
  - Initial implementation

Version 1.1
  - Implementation of subroutines
  - Better file IO
  - Cleanup of redundant code

Author: Jon Parsons
10-24-19

*/

#include <iostream>
#include <stdio.h>
#include "mole.h"

using namespace std;

int fileout(const vec& Arr, int m, int n, float dx, float dy, int f) {
/**
 Subroutine prints to file without replacement.
 Arguements: Array, array size in x, in y, step size in x, in y

*/

float x_cur, y_cur;

  FILE *outputFile = fopen("time_flow.dat","a");

  for(int x = 0; x < m; x++){
    x_cur = dx*((float) x);
    for(int y = 0; y < n; y++){
      y_cur = dy*((float) y);
      fprintf(outputFile, "%6.3f %6.3f %6.3f \n", x_cur, y_cur, Arr((x*n)+y));
    }
  }

 fprintf(outputFile, "X\n");

 fclose(outputFile);

 if (f == 1){
   FILE *frameFile = fopen("time_end.dat","a");

   for(int x = 0; x < m; x++){
     x_cur = dx*((float) x);
     for(int y = 0; y < n; y++){
       y_cur = dy*((float) y);
       fprintf(frameFile, "%6.3f %6.3f %6.3f \n", x_cur, y_cur, Arr((x*n)+y));
     }
   }

  fclose(frameFile);
 }

 return 0;

}

////////////////////////////////////////////////////////////////////////////////
int ini_fill(vec& Arr, int m, int n, float iCon){
/**
 Subroutine fills a vector with a constant initial value
 Arguements: Array, array size in x, in y, concentration value
*/

float bound = 5; // Determines when non-zero elements start away from y extremes

  // Input concentrations for t = 0
  for(int x = 0; x < m; x++){
    for(int y = 0; y < n; y++){
      if ((x == 0) && (y >= (bound-1)) && (y <= (n-bound))){
        Arr(x*(n)+y) = iCon;
      }
    }
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int main() {
  // Variables
  // Algorithm variables
  float t0, tf, dt; // initial, final, and spacing in time
  float dx, dy; // Change in space. Dimensions and grid chosen S.T. dx = dy
  float x0, xf, y0, yf; // Initial and final spatial coordinates
  int grdx, grdy, grdt; // Number of gridpoints in space and time
  int wrkgrdx, wrkgrdy, wrkgrdsq; // Working grids
  int ord = 2; // Order of accuracy
  // Equation Variables
  float Dc, Fv; // Dispersion coefficient, flow velocity
  float inC = 1.0; // Input concentration
  int flag; // Flag for printing individual frames

// Variable assignments
flag = 0;

t0 = 0.0;
tf = 0.2;
x0 = 0.0;
xf = 130.0;
y0 = 0.0;
yf = 32.5;

grdx = 30;
grdy = grdx; // square matrice
grdt = 200;

wrkgrdx = grdx+2;
wrkgrdy = grdy+2;
wrkgrdsq = wrkgrdx*wrkgrdy;

dx = (xf-x0)/(float)(grdx-1);
dy = (yf-y0)/(float)(grdy-1);
dt = (tf-t0)/(float)(grdt-1);

Dc = 75.0;
Fv = 15.0;

// Get mimetic operators
Gradient G(ord,grdx,grdy,dx,dy);
Divergence D(ord,grdx,grdy,dx,dy);
Interpol I(grdx,grdy,0.5,0.5);

// Create concentration vector
vec C(wrkgrdsq);
C.fill(0.0);

vec V(wrkgrdsq);
V.fill(Fv);

ini_fill(C,wrkgrdx,wrkgrdy,inC);

// Print initial values
fileout(C,wrkgrdx,wrkgrdy,dx,dy,flag);

// Iteration in time loop
for(int t = 1; t <= grdt; t++){

  if (t == grdt/2){
    flag = 1;
  }

  // differential operations
  C += dt*(D*(Dc*(G*C)) - D*(I*(V%C)));
  // output current values
  fileout(C,wrkgrdx,wrkgrdy,dx,dy,flag);

}

return 0;

}
