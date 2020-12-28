/**
Program written for the final project, COMP670.

This program will analyse a system of equations relating to a bacteria-phage
ecosystem using the Mimetic Operators Library Enhanced (MOLE)

Version 1.0
  - Initial implementation

Author: Jon Parsons
12-3-19

*/

#include <iostream>
#include <stdio.h>
#include "mole.h"

using namespace std;

// Commonly used variables
float Ad = 5.0E8; // inverse of typical value of 2.0E-8
float ra = 0.04; // hr^-1
float D_b = 0.01; // hr^-1

////////////////////////////////////////////////////////////////////////////////
int fileout(const vec& Arrp, vec& Arrb, int m, int n, float dx, float dy, int f) {
/**
 Subroutine prints to file without replacement.
 Arguements: Array, array size in x, in y, step size in x, in y

*/

float x_cur, y_cur;

  FILE *outputFile = fopen("time_flow_p.dat","a");

  for(int x = 0; x < m; x++){
    x_cur = dx*((float) x);
    for(int y = 0; y < n; y++){
      y_cur = dy*((float) y);
      fprintf(outputFile, "%6.3f %6.3f %6.3f \n", x_cur, y_cur, Arrp((x*n)+y));
    }
  }

 fprintf(outputFile, "X\n"); // Frame seperator

 fclose(outputFile);

 if (f == 1){
   FILE *frameFile = fopen("time_mid_p.dat","a");

   for(int x = 0; x < m; x++){
     x_cur = dx*((float) x);
     for(int y = 0; y < n; y++){
       y_cur = dy*((float) y);
       fprintf(frameFile, "%6.3f %6.3f %6.3f \n", x_cur, y_cur, Arrp((x*n)+y));
     }
   }

  fclose(frameFile);
 }

 FILE *outputFileb = fopen("time_flow_b.dat","a");

 for(int x = 0; x < m; x++){
   x_cur = dx*((float) x);
   for(int y = 0; y < n; y++){
     y_cur = dy*((float) y);
     fprintf(outputFileb, "%6.3f %6.3f %6.3f \n", x_cur, y_cur, Arrb((x*n)+y));
   }
 }

fprintf(outputFileb, "X\n"); // Frame seperator

fclose(outputFileb);

if (f == 1){
  FILE *frameFileb = fopen("time_mid_b.dat","a");

  for(int x = 0; x < m; x++){
    x_cur = dx*((float) x);
    for(int y = 0; y < n; y++){
      y_cur = dy*((float) y);
      fprintf(frameFileb, "%6.3f %6.3f %6.3f \n", x_cur, y_cur, Arrb((x*n)+y));
    }
  }

 fclose(frameFileb);
}

 return 0;

}

////////////////////////////////////////////////////////////////////////////////
int ph_fill(vec& Arr, vec cc, vec ba, int m, int n){
/**
 Subroutine fills the phage population vector
*/

for(int x = 0; x < m; x++){
  for(int y = 0; y < n; y++){
    Arr(x*n+y) = (ra*(1.0 - ba(x*n+y)/cc(x*n+y)) - D_b)*Ad;
  //  cout << (Arr(x*n+y)) << "  " << ba/cc(x*n+y) << "\n";
  }
}

return 0;

}

////////////////////////////////////////////////////////////////////////////////
int ba_fill(vec& Arr, vec dp, float bu, int m, int n){
/**
 Subroutine fills the phage population vector
*/

for(int x = 0; x < m; x++){
  for(int y = 0; y < n; y++){
    Arr(x*n+y) = dp(x*n+y)*Ad/bu;
  //  cout << (Arr(x*n+y)) << "  " << ba/cc(x*n+y) << "\n";
  }
}

return 0;

}

////////////////////////////////////////////////////////////////////////////////
float cc_fill(vec& Arr, float m, float n, float def, float ele){
/**
  Subroutine fills the initial distribution of carrying capacity
*/

for(int x = 0; x < m; x++){
  for(int y = 0; y < n; y++){
    if ((x > 10) && (x < 15) && (y > 10) && (y < 15)){
      Arr(x*n+y) = ele;
    } else if ((x > 75) && (x < 85) && (y > 75) && (y < 85)){
      Arr(x*n+y) = ele;
    } else if ((x > 35) && (x < 40) && (y > 25) && (y < 40)){
      Arr(x*n+y) = ele;
    } else if (Arr(x*n+y) < def) {
      Arr(x*n+y) = def;
    }
  }
}

return 0;

}

////////////////////////////////////////////////////////////////////////////////
float vel_fill(vec& Arr, float m, float n){
/**
  Subroutine to fill the velocity field vector
  Direction of flow in positive x and y
  This version decreases as y increases.
*/

for(int x = 0; x < m; x++){
  for(int y = 0; y < n; y++){
    Arr(x*n+y) =(n - y)*(0.3);
  }
}

}

////////////////////////////////////////////////////////////////////////////////
float dp_fill(vec& Arr, float m, float n, int f){
/**
  Subroutine to adjust the rate of phage death
*/

for(int x = 0; x < m; x++){
  for(int y = 0; y < n; y++){
    if (f == 1){
      if ((x > 12) && (x < 40)){
        Arr(x*n+y) = 0.01;
      } else {
        Arr(x*n+y) = 0.04;
      }
    } else {
      if (Arr(x*n+y) < 0.01){
        Arr(x*n+y) = 0.01;
      } else if (Arr(x*n+y) > 0.04){
        Arr(x*n+y) = 0.04;
      }
    }
  }
}

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
  float bu; // Burst size
  float D_b, D_p_std; // Natural death rate of bacteria, phage
  float amb, cor; // Sugar concentration for ambient and coral adjacent water
  float Dc; // Dispersion and velocity coefficients

  // State variables
  int flag; // Flag for printing individual frames

// Variable assignments
flag = 0;

t0 = 0.0;
tf = 0.3;
x0 = 0.0;
xf = 100.0;
y0 = 0.0;
yf = 100.0;

grdx = 100;
grdy = grdx; // square matrice
grdt = 300;

wrkgrdx = grdx+2;
wrkgrdy = grdy+2;
wrkgrdsq = wrkgrdx*wrkgrdy;

dx = (xf-x0)/(float)(grdx-1);
dy = (yf-y0)/(float)(grdy-1);
dt = (tf-t0)/(float)(grdt-1);

bu = 50.0;
D_p_std = 0.04;
amb = 0.42E7;
cor = 3.98E7;

Dc = 50;

// Get mimetic operators
Gradient G(ord,grdx,grdy,dx,dy);
Divergence D(ord,grdx,grdy,dx,dy);
Interpol I(grdx,grdy,0.5,0.5);

// Create concentration vectors
vec phage(wrkgrdsq);
phage.fill(0.0);
vec bact(wrkgrdsq);
bact.fill(0.0);
vec cc(wrkgrdsq);
cc.fill(amb);
vec vf(wrkgrdsq);
vf.fill(0.0);
vec Dp(wrkgrdsq);
Dp.fill(D_p_std);

// Initial filling of vectors
vel_fill(vf,wrkgrdx,wrkgrdy);
cc_fill(cc,wrkgrdx,wrkgrdy,amb,cor);
ph_fill(phage,cc,bact,wrkgrdx,wrkgrdy);
ba_fill(bact,Dp,bu,wrkgrdx,wrkgrdy);

// Print initial values
fileout(phage,bact,wrkgrdx,wrkgrdy,dx,dy,flag);

// Iteration in time loop
for(int t = 1; t <= grdt; t++){

  // Flag for single frame output
  if (t == grdt/2){
    flag = 1;
  } else {
    flag = 0;
  }

  // differential operations
  cc += dt*(D*(Dc*(G*cc)) - D*(I*(vf%cc)) - amb*vf);
  // Determine start of phage death change, before then, constant
  if (t == grdt/3){
  // Fun Fact: if you write t = n  instead of  t == n it assigns t to equal n
  // and makes an infinite loop
    dp_fill(Dp,wrkgrdx,wrkgrdy,1);
  }
 if (t >= grdt/3){
    Dp += dt*(D*(0.01*(G*Dp)) - D*(I*(vf%Dp)));
    dp_fill(Dp,wrkgrdx,wrkgrdy,0);
  }

  // output current values
  cc_fill(cc,wrkgrdx,wrkgrdy,amb,cor);

  // Data output
  fileout(phage,bact,wrkgrdx,wrkgrdy,dx,dy,flag);
  // Update phage and bacteria
  ph_fill(phage,cc,bact,wrkgrdx,wrkgrdy);
  ba_fill(bact,Dp,bu,wrkgrdx,wrkgrdy);
  // Refill sources
  cc_fill(cc,wrkgrdx,wrkgrdy,amb,cor);

}

return 0;

}
