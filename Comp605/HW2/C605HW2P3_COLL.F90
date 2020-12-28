! Code written for homework 2 problem 3, COMP605.
! Program calculates an integral using collective methods
! Program finds an approximation to pi using the integration
! of a function and compares to a known value.

! Author: Jon Parsons
! Date: 2-28-2019
module functions

contains

real (kind=8) function f(y)
! Function for use in integration

implicit none
	real (kind=8) :: y
	
	
	f = 4.0/(1.0 + y*y)

end function

end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program integralcoll

use mpi

implicit none
	real (kind=8)			:: a, b ! Global Endpoints
	real (kind=8)			:: local_a, local_b ! Local Endpoints
	real (kind=8)			:: h ! Stepsize
	integer					:: n ! Global number of subdivisions
	integer					:: local_n ! Local nuber of subdivisions
	integer					:: i ! Looping integer
	real (kind=8),parameter :: pi=acos(-1.0) ! A known, good value of pi
	real (kind=8)			:: results, local_results ! Global and local results of integration
	
	integer					:: world_size	! Holds number of processors available
	integer					:: mpi_err ! MPI error variable
	integer					:: local_num ! Holds local address of processor
	
	
	! Initialize MPI
	call MPI_INIT (mpi_err)
	call MPI_COMM_SIZE (MPI_COMM_WORLD,world_size, mpi_err)
	call MPI_COMM_RANK (MPI_COMM_WORLD,local_num, mpi_err)

! Endpoints
a = 0.0
b = 1.0

call input(local_num,world_size,a,b,n)

! Determine local subdivisions and stepsize
local_n = n/world_size

h = (b-a)/float(n)

! Initialize
results = 0.0
local_results = 0.0

! Determine local endpoints
local_a = a + local_num*float(local_n)*h
local_b = local_a + float(local_n)*h

call simpint(local_a,local_b,local_n,h,local_results)

! Collect results
call MPI_REDUCE(local_results,results,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)

! Output to user
if (local_num .eq. 0) then
	write(*,*) "Program Complete."
	write(*,*) "Calculated:", results
	write(*,*) "True:", pi
	write(*,*) "Delta:", results-pi
end if

call MPI_FINALIZE(mpi_err)
	
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine input(rank,size,left,right,num)
! Subroutine handles user input and data distribution to daughter processes

use mpi

	implicit none
		integer			:: rank, num, size
		real(kind=8)	:: left, right
		
		integer			:: i
		integer			:: ierr
		
if (rank .eq. 0) then
	write(*,*) "Please enter the exponent for use in finding number of sub-intervals (N=10^X):"
	read(*,*) num
	
	num = 10**num
end if

! MPI_BCAST sends data to all available processes
call MPI_BCAST(left,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(right,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(num,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

end subroutine
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
subroutine simpint(left,right,num,dx,out)
! Subroutine contains a Simpson's quadrature method that analyzes a function found in the 
! module 'functions'

use functions

implicit none
	real(kind=8)		:: left, right
	integer				:: num
	real(kind=8)		:: out, dx, xi
	
	integer				:: i
		
	out = f(left) + f(right)
	
	do i = 2, num-2, 2
		xi = left + float(i)*dx
		out = out + 2.0*f(xi) + 4.0*f(xi+dx)
	end do
	
	out = out*dx/3.0
	
end subroutine
		
	
		
		
		
		
		
		




