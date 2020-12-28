module params
! Holds variables from physical propoerties of 
! system

	implicit none
		real*8 	:: v, l
		real*8	:: w, tb
		real*8	:: td
		real*8	:: g, e
		real*8	:: eps0, epsr
		real*8	:: vcoeff, h
		real*8	:: dx
		

end module


program mems_sln
! Program written for final project, Comp521
! Program uses the Newton-Raphson iterative linear system
! solving method to find the deflection of a cantilever
! beam when differing force is applied in the form of 
! an electric field. 
! Author: Jon Parsons
! Date: 12-13-2018


use params

implicit none
! Variable vector, holding vector, error holding arrays
	real*8,allocatable	:: vars(:), varstmp(:), error(:,:)
! Hold the Jacobian and its inverse
	real*8,allocatable	:: jac(:,:), jacinv(:,:)
! Holds the evaluation of the system 
	real*8,allocatable	:: funcval(:)
! Number of divisions, current iteration, maximum iterations
	integer				:: grdpts, iter, maxiter
! Looping variables
	integer				:: i, j, k, n
! File Variables
	character*30		:: filename, errname
! Error handling variables
	real*8				:: errtol, errsum
! Step size in x
	real*8				:: xi

! Variable initializations
grdpts = 10
maxiter = 1000
errtol = 0.001
iter = 1

w = 1.0
tb = 3.0
l = 1.0
td = 10.0e3
g = 1.0
e = 0.169
eps0 = 8.85e-6
epsr = 11.68e6

vcoeff = ((eps0)/(2.0*e*((tb**3)/12.0)))
h = td/(g*epsr)

write(*,*) vcoeff, h

! File name format statements
51 format("Volts",1f5.0,"grid",1i2,".dat")
52 format("Err",1i2,"grd.dat")

! Outer loop iterates over step-size
do k = 1, 3, 1
	
	! Initialize variables that change with stepsize
	v = 0.0
	
	dx = l/float(grdpts)

	! Array allocation statements
	allocate(vars(grdpts+4))
	allocate(varstmp(grdpts+4))
	allocate(jac(grdpts+4,grdpts+4))
	allocate(jacinv(grdpts+4,grdpts+4))
	allocate(funcval(grdpts+4))
	allocate(error(grdpts+4,maxiter))
	
	! Loop iterates over voltages
	do j = 1, 7, 1
	
		! Initialize variables which change with voltage
		jacinv = 0.0
		jac = 0.0
		varstmp = 1.0
		vars = g
		
		! Iteration reset
		iter = 1

		! Store current variable vector (guess)
		error(:,iter) = vars
	
		! Initial evaluation at guess
		call funevals(vars,funcval,grdpts+4)
		
		! Determine the file names
		write(filename,51) v, grdpts
		write(errname,52) grdpts
		
		! Open the files
		open(unit=11,file=trim(filename),status="replace",position="append")
		open(unit=12,file=trim(errname),status="replace",position="append")
	
		! Main working loop, performs Newton-Raphson root-finding method
		do while ((maxval(abs(varstmp)).gt. errtol).and.(iter .lt. maxiter))
	
			call jacob(vars,jac,grdpts+4)
	
			call inverter(jac,jacinv,grdpts+4)	
			call matmulsq(jacinv,funcval,varstmp,grdpts+4)
			
			! Update variable vector
			vars = vars - varstmp
			
			! Evaluate for current variable vector
			call funevals(vars,funcval,grdpts+4)
		
			! Update iteration counter
			iter = iter + 1	
			
			! Store variable vector
			error(:,iter) = vars

		end do
		
		! Output to screen
		write(*,*) "Iterations:",iter, "V=",v,"grid:",grdpts
	
		! Update voltage
		v = v + 50.0
		
		! Output beam to file
		do i = 3, grdpts, 1
			xi = dx*float(i-3)
			write(11,*) 5.0*xi, 3.0*vars(i), funcval(i)
		end do
		
		! Calculate error and write to error file
		! file as : Iteration  Variable - root
		do i = 1, iter, 1
			
			errsum = 0.0
			
			do n = 1, grdpts+4, 1
			
				errsum = errsum + abs(error(n,i))
				
			end do

			write(12,*) i, log(abs(errsum - abs(sum(vars))))
			
		end do
		

		close(11)
		close(12)
		
	end do
	
	write(*,*)
	deallocate(vars)
	deallocate(varstmp)
	deallocate(jac)
	deallocate(jacinv)
	deallocate(funcval)
	deallocate(error)
		
	! Update gridpoints
	grdpts = grdpts + grdpts
	
end do

end program

subroutine funevals(in,out,dime)
! Function evaluates system of eqns for input vector of 
! variables

use params

implicit none
	integer		:: dime
	real*8		:: in(dime), out(dime)
	real*8		:: f1, f2, f3
	real*8		:: x, y, z
	integer		:: i
	
! Driving function
f1(x) = vcoeff*(v*v)*(dx**4)/((x+h)**2)

! Boundary conditions
out(1) = in(3) - 1.0
out(2) = in(1)-8.0*in(2)+8.0*in(4)-in(5)

! Evaluation in beam
do i = 3, dime-2, 1
	out(i) = in(i-2)-4.0*in(i-1)+6.0*in(i)-4.0*in(i+1)+in(i+2) + f1(in(i))
end do

! Boundary Conditions
out(dime-1) = -in(dime-4)+16.0*in(dime-3)-30.0*in(dime-2)+16.0*in(dime-1)-in(dime)
out(dime) = -in(dime-4)+2.0*in(dime-3)-2.0*in(dime-1)+in(dime)

end subroutine

subroutine jacob(in,out,dime)
! Evaluates the jacobian of the system of equations
! for input variable vector

	use params
	
	implicit none
		integer	:: dime
		real*8	:: in(dime), out(dime,dime)
		integer :: i
		real*8	:: f, x
	
! Driving function 
f(x) = -vcoeff*(v*v)*2.0*(dx**4)/((x+h)**3)
		
! Initialize
out = 0.0

! Boundary conditions
out(1,3) = 1.0
out(2,1) = 1.0
out(2,2) = -8.0
out(2,4) = 8.0
out(2,5) = -1.0

! Working loop for beam values
do i = 3, dime-2 ,1
		out(i,i-2) = 1.0
		out(i,i-1) = -4.0
		out(i,i) = 6.0 + f(in(i))
		out(i,i+1) = -4.0
		out(i,i+2) = 1.0
end do

! Boundary conditions
out(dime-1,dime-4) = -1.0
out(dime-1,dime-3) = 16.0
out(dime-1,dime-2) = -30.0
out(dime-1,dime-1) = 16.0
out(dime-1,dime) = -1.0

out(dime,dime-4) = -1.0
out(dime,dime-3) = 2.0
out(dime,dime-1) = -2.0
out(dime,dime) = 1.0
		
end subroutine		
		

subroutine inverter(a,ainv,dime)
! Function inverts input matrix a to ainv
! Partial pivoting applied

implicit none
	integer				:: dime
	real*8				:: a(dime,dime), ainv(dime,dime)
	integer				:: i,j, dub, p
	real*8,allocatable	:: work(:,:)
	real*8				:: t

dub = 2*dime

a(1,1) = 0.00000001
ainv = 0.0


allocate(work(dime,dub))

do j = 1, dime, 1
	do i = 1, dime, 1
		work(i,j) = a(i,j)
		work(i,j+dime) = 0.0
		if (i .eq. j) then
			work(i,j+dime) = 1.0
		end if
	end do
end do

do i = 1, dime, 1

	p = i
		
	do j = i+1, dime, 1
		if (abs(a(j,i)).gt.abs(a(p,i))) then
			p = j
		end if
	end do
	
	if (p .ne. i) then
		do j = 1, dime, 1
			t = a(i,j)
			a(i,j) = a(p,j)
			a(p,j) = t
		end do
	end if
	
	work(i,i+1:dub) = work(i,i+1:dub)/work(i,i)
	work(i,i) = 1.0
	do j = 1, dime, 1
		if (j .eq. i) then
			cycle
		end if
	
		work(j,i+1:dub) = work(j,i+1:dub)-work(i,i+1:dub)*work(j,i)
		work(j,i) = 0.0
	end do
end do

do i = 1, dime, 1
	ainv(:,i) = work(:,i+dime)
end do

deallocate(work)

end subroutine


subroutine matmulsq(insq,inv,out,dime)
! Multiplies a square matrix and a vector
	implicit none
		integer		:: dime
		real*8		:: insq(dime,dime),inv(dime)
		real*8		:: out(dime)
		integer		:: i, j
		
	out = 0.0
	
do i = 1, dime, 1
	
	do j = 1, dime, 1
		
		out(i) = out(i) + insq(i,j)*inv(j)
		
	end do
	
end do

end subroutine
		





















