! Code written for homework 2 problem 3, COMP605.
! Program calculates an integral using point-to-point mpi protocols
! The integral finds an approximation to pi and compares to a known 
! close value.

! Author: Jon Parsons
! Date: 2-28-2019

module functions

contains

real (kind=8) function f(y)
! Function contains the formula leading to the approximation of pi

implicit none
	real (kind=8) :: y
	
	
	f = 4.0/(1.0 + y*y)

end function

end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program integralp2p

use mpi

implicit none
	real (kind=8)			:: a, b ! Global endpoints
	real (kind=8)			:: local_a, local_b ! Local endpoints
	real (kind=8)			:: h ! Stepsize
	integer					:: n ! Global number of subdivisions
	integer					:: local_n ! Local number of subdivisions
	integer					:: i ! Looping integer
	real (kind=8),parameter :: pi=acos(-1.0) ! This method gives a good approximation of pi to check against
	real (kind=8)			:: results, local_results ! Global and local values of integration
	
	integer					:: world_size	! Holds number of processors being used
	integer					:: mpi_err ! Error checking variable
	integer					:: local_num ! Holds local address of processor
	
	! Initialize MPI
	call MPI_INIT (mpi_err)
	call MPI_COMM_SIZE (MPI_COMM_WORLD,world_size, mpi_err)
	call MPI_COMM_RANK (MPI_COMM_WORLD,local_num, mpi_err)

! Initialize endpoints
a = 0.0
b = 1.0

call input(local_num,world_size,a,b,n)

! Determine local number of subdivisions
local_n = n/world_size

! Determine stepsize
h = (b-a)/float(n)

results = 0.0

local_results = 0.0

! Determine the endpoints each processor will be using
local_a = a + local_num*float(local_n)*h
local_b = local_a + float(local_n)*h

call simpint(local_a,local_b,local_n,h,local_results)

! Returns the values each local calculation found to the master for output.
if (local_num .ne. 0) then
	call MPI_SEND(local_results,1,MPI_REAL8,0,0,MPI_COMM_WORLD,mpi_err)
else
	results = results + local_results
	do i = 1, world_size-1, 1
	call MPI_RECV(local_results,1,MPI_REAL8,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_err)
		results = results + local_results
	end do
end if 

! Return final results to user.
if (local_num .eq. 0) then
	write(*,*) "Program Complete."
	write(*,*) "Calculated:", results
	write(*,*) "True:", pi
	write(*,*) "Delta:", results-pi
end if

call MPI_FINALIZE(mpi_err)
	
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine input(rank,size,left,right,num)
! This subroutine streamlines the user input process and handles 
! distribution to daughter processes.

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
	 
	! Distributes values to each processor individually
	do i = 1, size-1, 1
		call MPI_SEND(right,1,MPI_REAL8,i,0,MPI_COMM_WORLD,ierr)
		call MPI_SEND(left,1,MPI_REAL8,i,0,MPI_COMM_WORLD,ierr)
		call MPI_SEND(num,1,MPI_INT,i,0,MPI_COMM_WORLD,ierr)		
	end do
	
else 
	! Each daughter process receives the data here.
	call MPI_RECV(right,1, MPI_REAL8,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	call MPI_RECV(left,1,MPI_REAL8,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
	call MPI_RECV(num,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
subroutine simpint(left,right,num,dx,out)
! This subroutine utilizes the Simpson's quadrature method to 
! calculate the integral of a function found in module 'functions'

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
		
	
		
		
		
		
		
		




