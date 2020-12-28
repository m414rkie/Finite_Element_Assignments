/**
Program written for Homework 7, MATH693B.

This program will solve the problem 10.4.1 from Strikwerda using mimetic
operators. The problem involves the solution of the heat equation with two sets
of initial conditions and periodic boundaries. The nature of mimetic operators
will prevent the periodic boundaries, but no accuracy will be lost. The program
will implement the Forward-Time scheme for evolution in time with the step-size
k == h^2 in order to preserve second order accuracy.

Author: Jon Parsons
5-7-2020

*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include "mole.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int fileout(const vec& Arre, const vec& Arrn, int m, float dx, int flg, int hc)
{
/**
 Subroutine prints to file without replacement.
 Arguements:
 Arre - Array with exact wave values
 Arrn - Array with numerical wave values
 m    - Size of array
 dx   - step size
 flg  - determines which file to print to
 hc   - which iteration of h

*/

// Time series output
std::string fname1 = std::string("IC1_") + std::to_string(hc) + ".dat";
std::string fname2 = std::string("IC2_") + std::to_string(hc) + ".dat";

const char *fstring1 = fname1.c_str();
const char *fstring2 = fname2.c_str();

cout.precision(5);

if (flg == 1){
  ofstream fileo;
  fileo.open(fstring1, std::ios_base::app);

for(int x = 0; x < m; x++){
  float x_cur = -1.0 + dx*((float) x);
   fileo << x_cur << " " << Arre(x) << " " << Arrn(x) << endl;
}

fileo << "X\n"; // Frame seperator

} else {
  ofstream fileo;
  fileo.open(fstring2, std::ios_base::app);

for(int x = 0; x < m; x++){
  float x_cur = -1.0 + dx*((float) x);
   fileo << x_cur << " " << Arre(x) << " " << Arrn(x) << endl;
}

fileo << "X\n"; // Frame seperator

}

return 0;

}

////////////////////////////////////////////////////////////////////////////////
int ini_fill_a(vec& Arr, int m, float dx){
/**
 Subroutine fills the initial values for the first set of initial conditions
  Arguements:
  Arr - Array with wave values
  m   - Size of array
  dx  - step size

*/

for(int x = 0; x < m; x++){
   float x_c = -1.0 + dx*((float) x);
    if (abs (x_c) < 0.5) {
        Arr(x) = 1.0;
    } else if (abs (x_c) == 0.5){
        Arr(x) = 0.5;
    } else {
    Arr(x) = 0.0;
  }
}

return 0;

}

////////////////////////////////////////////////////////////////////////////////
int ini_fill_b(vec& Arr, int m, float dx){
/**
 Subroutine fills the initial values for the second set of initial conditions
  Arguements:
  Arr - Array with wave values
  m   - Size of array
  dx  - step size

*/

float pi = 3.14159265;

for (int x = 0; x < m; x++){
  float x_c = -1 + x*dx;
  Arr(x) = cos (pi*x_c);
}

return 0;

}

////////////////////////////////////////////////////////////////////////////////
float exact_a(float x, float t){
/**
Subroutine containing exact solution of the first set of initial conditions
Arguements:
 Arr - Array with wave values
 m   - Size of array
 dx  - step size
 t   - time

*/

// summation limit
int l = 50;

float pi = 3.14159265;

float sum = 0;
for (int k = 0; k <= l; k++){
  float n_one = pow(-1,k);
  float ll = float (2*k) + 1;
  sum += exp (-t*ll*ll*pi*pi) * (n_one/ll) * cos (ll*pi*x);
}

float sln = (2/pi)*sum + 0.5;

return sln;

}

////////////////////////////////////////////////////////////////////////////////
float exact_b(float x, float t){
/**
Subroutine containing exact solution of the second set of initial conditions
Arguements:
 Arr - Array with wave values
 m   - Size of array
 dx  - step size
 t   - time

*/

float pi = 3.14159265;

float sln = cos (pi*x) * exp (-pi*pi*t);


return sln;

}

////////////////////////////////////////////////////////////////////////////////
float L2_norm(vec Num, vec Exa, float h){
/**
Subroutine to define the L2 norm
 Num - values of the numerically calculated solution
 Exa - values of the exact solutions
 h   - stepsize

*/

int n = Num.size();

float sum = 0;
for (int i = 0; i < n; i++){
  sum += (Exa(i) - Num(i))*(Exa(i) - Num(i));
}

float L2_n = sqrt (sum*h);

return L2_n;

}

////////////////////////////////////////////////////////////////////////////////
float acc_ord(vec nrms){
/**
Subroutine to determine the order of accuracyas determined by the norms
 nrms - vector containing the calculated norms
*/

float sum = 0;
int n = nrms.size();
float nn = 0.0;

for (int i = 0; i < n-1; i++){
  sum += log2 (nrms(i)/nrms(i+1));
  nn += 1.0;
}

float ord = sum/nn;

return ord;

}

////////////////////////////////////////////////////////////////////////////////
int main() {
  // Variables
  // Algorithm variables
  float ti, tf, dt; // initial, final, and spacing in time
  float mu; // dt/dx^2
  float dx; // Change in space. Dimensions and grid chosen S.T. dx = dy
  float xi, xf; // Initial and final spatial coordinates
  int grdx, grdt; // Number of gridpoints in space and time
  int wrkgrdx; // Working grids
  int ord = 2; // Order of accuracy


  // State variables
  int flag; // Flag for printing individual frames

// Variable assignments
mu = 0.4;
ti = 0.0;
tf = 1.0;
xi = -1.0;
xf = 1.0;

dx = 1.0/10.0;

vec l2_1(3);
l2_1.fill(0.0);
vec l2_2(3);
l2_2.fill(0.0);

for (int h = 0; h <= 2; h++){
  cout << "Starting an h " << dx << endl;
  dt = dx*dx*mu;

  grdx = round ((xf - xi)/dx);
  grdt = round ((tf - ti)/dt);

  // Get mimetic operators
  Gradient G(ord,grdx,dx);
  Divergence D(ord,grdx,dx);


  // Create wave vectors
  vec I1_a(grdx+2);
  I1_a.fill(0.0);
  vec I2_a(grdx+2);
  I2_a.fill(0.0);
  vec I1_u(grdx+2);
  I1_u.fill(0.0);
  vec I2_u(grdx+2);
  I2_u.fill(0.0);

  // Initial filling of vectors
  ini_fill_a(I1_a,grdx+2,dx);
  ini_fill_a(I1_u,grdx+2,dx);
  ini_fill_b(I2_a,grdx+2,dx);
  ini_fill_b(I2_u,grdx+2,dx);

  // Print initial values
  fileout(I1_a,I1_u,grdx+2,dx,1,h);
  fileout(I2_a,I2_u,grdx+2,dx,2,h);

  // Iteration in time loop
  for(int t = 1; t <= grdt; t++){
    float t_c = ti + float (t) * dt;
    // differential operations
      I1_u += dt*(D*G*(I1_u));
      I2_u += dt*(D*G*(I2_u));

    // Boundaries
    I1_u(0) = exact_a(xi,t_c);
    I2_u(0) = exact_b(xi,t_c);
    I1_u(grdx+1) = exact_a(xf,t_c);
    I2_u(grdx+1) = exact_b(xf,t_c);

    // exact updates
    for (int x = 0; x < grdx+2; x++){
      float x_c = xi + dx* float (x);
      I1_a(x) =  exact_a(x_c,t_c);
      I2_a(x) =  exact_b(x_c,t_c);
    }
    if (t%10 == 0){
      fileout(I1_a,I1_u,grdx+2,dx,1,h);
      fileout(I2_a,I2_u,grdx+2,dx,2,h);
    }
  }
  cout << "finding L2's " << h << endl;
  l2_1(h) = L2_norm(I1_u,I1_a,dx);
  l2_2(h) = L2_norm(I2_u,I2_a,dx);

  dx /= 2;
}

float i1_order = acc_ord(l2_1);
float i2_order = acc_ord(l2_2);

cout << "L2 Norms and Order of Accuracy\n";
cout << "First set of Initial Conditions\n";
cout << l2_1 << endl;
cout << i1_order << endl;
cout << "Second set of Initial Conditions\n";
cout << l2_2 << endl;
cout << i2_order << endl;

return 0;

}
