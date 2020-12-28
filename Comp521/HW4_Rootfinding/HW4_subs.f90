!Subroutines and functions for HW 4 in Comp 521.
! Author: Jon Parsons
! Date: 11-1-2018

module functions

contains

real function fun(x)
! Function containing function to be evaluated.
	real	:: pi = acos(-1.0)
	real	:: x
	
fun = exp(-0.2*x)*sin(0.2*pi*x)

end function fun

real function df(x)
! Function containing derivative of function of interest
	real	:: pi = acos(-1.0)
	real	:: x
	
df = -0.2*exp(-0.2*x)*(sin(0.2*pi*x) - pi*cos(0.2*pi*x))

end function df

end module

subroutine bisec(a,b,errtol)
! Subroutine containing code for bisection method of rootfinding.
use functions

implicit none
integer				:: i, k	! Looping integers
real				:: a, b, errtol, currval ! Range values, tolerance, and current guess
integer				:: maxiter ! Max number of iterations 
real				:: tryval ! Holding variable for current value
real, allocatable 	:: holdarr(:,:) ! Array for holding values at each iteration
real				:: rootval ! Final value

! Initializations and allocations
maxiter = floor(sqrt((b-a)/errtol))
allocate(holdarr(3,maxiter))
holdarr = 0.0

! Output file. Column 1: Iteration number
! Column 2: log value of difference between current value and found value
! Column 3: log value of range size
open(unit=15,file="Bisectvals.dat",status="replace",position="append")

do i = 1, maxiter, 1
	
	tryval = (b+a)/2.0
	currval = fun(tryval)
	holdarr(1,i) = float(i)
	holdarr(2,i) = tryval
	holdarr(3,i) = (b-a)
	
	! Logic statements. If current value does not result in an acceptable value the range is shortened
	! and the current value is updated.
	if ((abs(currval) .le. errtol).or.(((b-a)/2.0) .lt. errtol)) then
		rootval = tryval
		goto 101
	else if (((currval.gt.0.0).and.(fun(a).gt.0.0)).or.((currval.lt.0.0).and.(fun(a).lt.0.0))) then
		a = tryval
	else
		b = tryval
	end if
		
end do
	
if (i .eq. maxiter) then
	write(*,*) "Bisection method failed to converge"
	exit
end if
	
101	write(*,*) "Bisection method result: ", rootval

do k = 1, i, 1
	write(15,*) holdarr(1,k), log(abs(holdarr(2,k)-rootval)), log(holdarr(3,k))
end do

close(15)
deallocate(holdarr)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine regfals(a,b,errtol)
! Contains code for the false position method.

use functions

implicit none
	real				:: a, b, errtol ! Range and tolerance
	integer				:: i, k, maxiter ! Looping integers and max. number of iterations
	real				:: citer, x, y, z, w ! Function variables
	real				:: rootval	! Final value
	real, allocatable	:: holdarr(:,:) ! Holds values at each iteration

! Function for updating current value.
citer(x,y,z,w) = y - (z*(y-x))/(z-w)

! Initializations and allocations.
maxiter = floor(sqrt((b-a)/errtol))
allocate(holdarr(3,maxiter))
holdarr = 0.0

! File for output.
! Column 1 : Iteration number
! Column 2 : Log of difference between current value and found value.
! Column 3 : Log of range size.
open(unit=15,file="Regfalsivals.dat",status="replace",position="append")

do i = 1, maxiter, 1
	
	rootval = citer(a,b,fun(b),fun(a))
	
	holdarr(1,i) = i
	holdarr(2,i) = rootval
	holdarr(3,i) = (b-a)
	
	if (abs(fun(rootval)) .le. errtol) then
		goto 102
	else if (((rootval.gt.0.0).and.(fun(a).gt.0.0)).or.((rootval.lt.0.0).and.(fun(a).lt.0.0))) then
		b = rootval
	else
		a = rootval
	end if
	
end do

if (i .eq. maxiter) then
	write(*,*) "False position method failed to converge"
	exit
end if

102	write(*,*) "Regula Falsi method result: ", rootval

do k = 1, i, 1
	write(15,*) holdarr(1,k), log(abs(holdarr(2,k) - rootval)), log(holdarr(3,k))
end do

close(15)
deallocate(holdarr)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine newtraph(a,b,errtol)
! Contains code for Newton-Raphson method of rootfinding

use functions

implicit none
	real				:: a, b, errtol ! Range and tolerance
	integer				:: i, k, maxiter ! Looping integers and max. number of iterations
	real, allocatable	:: holdarr(:,:) ! Holds values from each iteration
	real				:: rootval ! Final value
	
! Initializations and allocations
maxiter = floor(sqrt((b-a)/errtol))
allocate(holdarr(2,maxiter))
holdarr = 0.0

! Files for output
! Unit 15
! Column 1 : Iteration number
! Column 2 : log of difference between current value and root value
! Unit 16 contains updating function of this method
open(unit=15,file="NewtRaphsvals.dat",status="replace",position="append")
open(unit=16,file="NewtRaphsupfcn.dat",status="replace",position="append")

rootval = 2.75

do i = 1, maxiter, 1
	
	holdarr(1,i) = i
	holdarr(2,i) = rootval
	
	if (abs(fun(rootval)) .le. errtol) then
		goto 103
	else
		rootval = rootval - fun(rootval)/df(rootval)
	end if
	
end do

if (i .eq. maxier) then
	write(*,*) "Newton-Raphson method failed to converge."
	exit
end if

103 write(*,*) "Newton-Raphson result: ", rootval

do k = 1, i, 1
	write(15,*) holdarr(1,k), log(abs(holdarr(2,k) - rootval))
	write(16,*) log(holdarr(2,k)), log(holdarr(2,k+1))
end do

close(15)
close(16)
deallocate(holdarr)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine secmeth(a,b,errtol)
! Contains code for secant method of rootfinding.

use functions

implicit none
	real				:: a, b, errtol ! Range and tolerance
	integer				:: i, k, maxiter ! Looping integers and max. number of iterations
	real, allocatable	:: holdarr(:,:) ! Holds values at each iteration
	real				:: val3, val1, val2 ! Val1 and val2 are previous guess, val3 is current guess.
	
! Initial guesses
val1 = 2.75
val2 = 2.85

maxiter = floor(sqrt((b-a)/errtol))
allocate(holdarr(2,maxiter))
holdarr = 0.0

! Files for output. 
! Unit 15
! Column 1 : Iteration number
! Column 2 : log of difference between current value and found value
! Unit 16 is updating function of this method.
open(unit=15,file="secantvals.dat",status="replace",position="append")
open(unit=16,file="secantupfcn.dat",status="replace",position="append")

do i = 1, maxiter, 1
	
	val3 = val2 - fun(val2)*((val2-val1)/(fun(val2)-fun(val1)))
	
	holdarr(1,i) = i
	holdarr(2,i) = val3
	
	if (abs(fun(val3)) .le. errtol) then
		goto 104
	else
		val1 = val2
		val2 = val3
	end if
	
end do
	
if (i .eq. maxiter) then
	write(*,*) "Secant method failed to converge."
	exit
end if
	
104 write(*,*) "Secant Method result: ", val3

do k = 1, i, 1
	write(15,*) holdarr(1,k), log(abs(holdarr(2,k) - val3))
	write(16,*) log(holdarr(2,k)),log( holdarr(2,k+1))
end do

close(15)
deallocate(holdarr)

end subroutine
