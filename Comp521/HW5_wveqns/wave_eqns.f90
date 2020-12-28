program waves
! Program written for assignment 6 in COMP 521.
! Program solves three wave equations, each with
! a unique method. 



call prob1
call prob2
call prob3

end program

subroutine prob1
! Subroutine solves a wave equation using the Lax-Friedrichs method.

implicit none 
! Holds wave at all times (x,t), holds h and error metric.
	real, allocatable 	:: sln(:,:), lnerrvals(:,:)
! Discretization variables
	real				:: h, dt, xi, ti
! Domain variables
	real				:: ax, bx, at ,bt
! Number of grid points in x and time
	integer				:: grdptsx, grdptst
! Function variabls
	real				:: f, x, y, g
! Error variables
	real				:: slope, diff
! Looping integers
	integer				:: i, j ,k
	
! Wave at t=0
f(x) = 1.0 - abs(x)

! Initialize variables
ax = -2.0
bx = 10.0

at = 0.0
bt = 8.0

h = 0.1

dt = 0.001

! Allocate error array
allocate(lnerrvals(2,5))

! Open error file
open(unit=10,file="lnerr1.dat",status="replace",position="append")

do i = 1, 5, 1

	! Determine num. of steps at current h
	grdptsx = nint((bx-ax)/h)
	grdptst = nint((bt-at)/dt)
	
	! Allocate solution array
	allocate(sln(grdptsx,grdptst))
	
	sln = 0.0
	
	! Set t=0 for first point
	if (abs(ax) .lt. 1.0) then
		sln(1,1) = f(ax)
	end if
	
	! Set rest of t=0
	do j = 2, grdptsx, 1
		
		xi = ax + float(j)*h
		
		if (abs(xi) .lt. 1.0) then
			sln(j,1) = f(xi)
		end if
		
	end do
	
	! Set x=ax, x=bx
	do j = 2, grdptst, 1
		
		ti = at + float(j)*dt
		
		if (abs(ax-ti) .lt. 1.0) then
			sln(1,j) = f(ax-ti)
		end if
		
		if (abs(bx-ti) .lt. 1.0) then
			sln(grdptsx,j) = f(bx-ti)
		end if
		
	end do
	
	! Update through the rest of time
	do j = 2, grdptst, 1
		
		do k = 2, grdptsx-1, 1
			
			sln(k,j) = 0.5*(sln(k-1,j-1)+sln(k+1,j-1))-0.5*(dt/h)*(sln(k+1,j-1)-sln(k-1,j-1))
		
		end do
	
	end do
	
	diff = 0.0
	
	! Find differences
	diff = abs(f(ax-bt) - sln(1,grdptst))
	
	do j = 2, grdptsx, 1
		
		xi = ax + float(j)*h
		
		diff = diff + abs(f(xi-bt) - sln(j,grdptst))
		
	end do
		
	! Print error stuff
	lnerrvals(1,i) = log(h)
	lnerrvals(2,i) = log(diff)
	write(10,*) lnerrvals(1,i), lnerrvals(2,i)

	! Output wave at final step
	if (i .eq. 5) then
		
		open(unit=11,file="waves1.dat",status="replace",position="append")
		write(11,*) "x num. exact"
		
		do j = 1, grdptsx, 1
			
			xi = ax + float(j)*h
		
		if (abs(xi-bt) .lt. 1.0) then	
			write(11,*) xi, sln(j,grdptst), f(xi-bt)
		else
			write(11,*) xi, sln(j,grdptst), 0.0
		end if
		
		end do
		
		slope = (lnerrvals(2,5)-lnerrvals(2,1))/(lnerrvals(1,5)-lnerrvals(1,1))
		
		write(*,*) "Error slope of problem 1: ", slope
	
		close(11)
		
	end if
	
	! Update h
	h = h*0.5
	
	deallocate(sln)
	
end do

close (10)

end subroutine
	
subroutine prob2
! Solves a wave equation using explicit method


implicit none 
! Holds wave at each timestep (x,t), holds error values
	real, allocatable	:: sln(:,:), lnerrvals(:,:)
! current x and boundary values
	real				:: xi, ax, bx, at, bt
! step size variables
	real				:: h, k, r
! Error variables
	real				:: slope, diff
! pi
	real				:: pi = acos(-1.0)
! function variables
	real				:: f, g, x, y
! Number of gridpoints in time and space
	integer				:: grdptst, grdptsx
! Looping integers
	integer				:: i, j, l

! f(x) initial values, g(x,y) exact solution
f(x) = sin(pi*x) + sin(3.0*pi*x)
g(x,y) = sin(pi*x)*exp(-pi*pi*y) + sin(3.0*pi*x)*exp(-9.0*pi*pi*y)

! Value initialization
ax = 0.0
bx = 1.0
	
at = 0.0
bt = 0.1

h = 0.2
k = 0.02
	
r = k/(h*h)
	
! Prepare error file and array
open(unit=10,file="lnerr2.dat",status="replace",position="append")

allocate(lnerrvals(2,5))

do i = 1, 5, 1

	! Determine number of gridpoints 
	grdptsx = nint((bx-ax)/h)
	grdptst = nint((bt-at)/k)
	
	! Allocation of solution array
	allocate(sln(grdptsx,grdptst))
	
	sln = 0.0
	
	! Initialize t=0
	do j = 2, grdptsx-1, 1
	
		xi = ax + float(j)*h
		
		sln(j,1) = f(xi)
	
	end do
	
	! Update through time
	do j = 2, grdptst, 1
		
		do l = 2, grdptsx-1, 1
			
			sln(l,j) = sln(l,j-1) + r*(sln(l-1,j-1) - 2.0*sln(l,j-1) + sln(l+1,j-1))
			
		end do
		
	end do
	
	diff = 0.0
	
	! Determine errors
	do j = 1, grdptsx-1, 1
		
		xi = ax + float(j)*h
	
		diff = diff + (abs(g(xi,bt) - sln(j,grdptst)))
		
	end do
	
	! Output errors
	lnerrvals(1,i) = log(h)
	lnerrvals(2,i) = log(diff)
	write(10,*) lnerrvals(1,i), lnerrvals(2,i)
			
	! File output at final h value
	if (i .eq. 5) then
		
		open(unit=11,file="waves2.dat",status="replace",position="append")
		write(11,*) "x num. exact"
		
		do j = 1, grdptsx, 1
			
			xi = ax + float(j)*h

			write(11,*) xi, sln(j,grdptst), g(xi,bt)
		
		end do
		
		slope = (lnerrvals(2,5)-lnerrvals(2,1))/(lnerrvals(1,5)-lnerrvals(1,1))
		
		write(*,*) "Error slope of problem 2: ", slope
	
		close(11)
		
	end if
		
	! Update values
	h = h*0.5
	k = 0.5*h*h
	r = k/(h*h)
	
	deallocate(sln)
	
end do

close (10)
	
end subroutine
	
subroutine prob3
! Solves a wave equation using a successive relaxation method

implicit none
! Main array, working array, difference between working and main, error array.
	real, allocatable	:: sln(:,:), slnn(:,:), errar(:,:), lnerrvals(:,:)
! Boundary values
	real				:: ax, bx, ay, by
! Step sizes, relaxation parameter
	real				:: h, rho
! Error variables
	real				:: err, errtol
! Maximum iterations
	integer				:: numit
! Number of gridpoints
	integer				:: grdptsx, grdptsy
! Looping integers
	integer				:: i, j, k, c

! Initializations
ax = 0.0
bx = 1.0

ay = 0.0 
by = 1.0

h = 0.1

errtol = 0.000001

numit = nint(h/errtol)

! Allocate error array
allocate(lnerrvals(2,numit))


do i = 1, 2, 1	
	
	! Determine number of gridpoints
	grdptsx = nint((bx-ax)/h)
	grdptsy = nint((by-ay)/h)
	
	rho = 1.0

	! Allocate arrays
	allocate(sln(grdptsx,grdptsy))
	allocate(slnn(grdptsx,grdptsy))
	allocate(errar(grdptsx,grdptsy))

	lnerrvals = 0.0

	! Initialize arrays 
	sln = (10.0+50.0+100.0)*0.25
	slnn = 0.0
	! Boundary values
	sln(:,1) = 10.0
	sln(1,:) = 50.0
	sln(:,grdptsy) = 100.0
	sln(grdptsx,:) = 0.0
	
	! Iteration loop
	do c = 1,  numit, 1
		
		! set holding array to main array
		slnn = sln
	
		! Iteration through y, x
		do j = 2, grdptsy-1, 1
		
			do k = 2, grdptsx-1, 1
		
				sln(k,j) = 0.25*rho*(sln(k+1,j) + sln(k-1,j) + sln(k,j+1) + sln(k,j-1))

		end do
		
		end do
		
		! Error finding
		errar = sln - slnn
		
		err = sum(abs(errar))
	
		lnerrvals(1,c) = c
		lnerrvals(2,c) = abs(err)
	
		! Determing number of iterations
		if ((abs(err)) .le. (errtol)) then
			write(*,*) "h =",h,"Number of iterations: ",c		
			exit
		end if
		
		slnn = 0.0
		
		! Failure logic
		if (c .eq. numit) then
			write(*,*) "Prob. 3 did not converge at h = ",h
		end if
		
	end do

	! Update h
	h = h*0.1

	if (i .ne. 2) then
		deallocate(sln)
		deallocate(slnn)
		deallocate(errar)
	end if
	
end do

! File output
open(unit=10,file="errvals3.dat",status="replace",position="append")
	write(10,*) "ln(h)  ln(diff)"
	
	do j = 1, c, 1
		write(10,*) lnerrvals(1,j), lnerrvals(2,j)
	end do
	
	close(10)

open(unit=11,file="SORfin.dat",status="replace",position="append")

	do j = 1, grdptsy, 1
	
		do k = 1, grdptsx, 1
			
			write(11,*) j, k, sln(j,k)
			
		end do

	end do
	
close(11)
	
end subroutine




	
	
	
	
	
	
	
	
	
	
	