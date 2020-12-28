! Program for assignment 4 in Comp. 521
! Program written to compare root-finding methods. Calls functions which contain
! code suitable for testing the bisection method, the false position method, the 
! Newton-Raphson method, and the Secant method.
! 
! Author: Jon Parsons
! 11-1-2018
!
! Compilation notes: Compiles with HW4_subs.f90

program roots

use Rootfinding

implicit none
real		:: a, b
real		:: tol


! Bisection section

a = 2.5
b = 6.4
tol = 10e-6

call bisec(a,b,tol)

a = 2.5
b = 6.4
tol = 10e-6

call regfals(a,b,tol)

a = 2.5
b = 6.4
tol = 10e-6

call newtraph(a,b,tol)

a = 2.5
b = 6.4
tol = 10e-6

call secmeth(a,b,tol)


end program
