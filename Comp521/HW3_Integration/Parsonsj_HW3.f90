program difference

!	Written for use in completing homework three for Comp521.
! Program calculates the integral of a function using the 
! trapezoidal rule as well as the Simpsons' 1/3 rule. The 
! integration will be done at varying step sizes. Absolute
! and relative errors will be calculated at each step size 
! for each method and output to screen. A file will also 
! be produced containing the log-log values of both error
! and step size for plotting.
!
! Author: Jon Parsons, 10-16-2018

implicit none
! Common variables
	integer				:: i, j		 		 ! Looping integers
	real				:: f, x		 		 ! Function variables
	real				:: era, err		 	 ! Error function variables
	real, parameter		:: pi = acos(-1.0) 	 ! Pi
! Array Variables
	integer				:: arrsize			 ! Holds size of array at each step
	real, allocatable	:: valt(:),vals(:)	 ! Holds the value of the function at each step
! Integration Variables	
	real				:: dx, xi			 ! Step size and current x
	real				:: a, b				 ! Bounds
	real				:: intval	 		 ! Holds final result
	real, parameter		:: truval = 0.1446445197	! True value of integral
	
! Functions for evaluation and error finding
f(x) = exp(-2.0*x)*sin(2.0*pi*x)
era(x) = abs(x - truval)
err(x) = abs(x - truval)/truval

! Bounds of integration
a = 0.0
b = 3.5

! Format statements for output
51 format("Abs. error for ", i3," segments in Trap. :", f13.11)
52 format("Rel. error for ", i3," segments in Trap. :", f13.11)
53 format("Abs. error for ", i3," segments in Simp. :", f13.11)
54 format("Rel. error for ", i3," segments in Simp. :", f13.11)

! File opening for error plotting
open(unit=11,file="Errplottrap.dat",status="replace",position="append")
open(unit=12,file="Errplotsimp.dat",status="replace",position="append")

! Initial allocation statements
allocate(valt(2))
allocate(vals(3))

! Initial (single-segment) values input by hand for simplicity.
valt(1) = f(a)
valt(2) = f(b)

call trapint(3.5,2,vals,intval)

write(*,*) "Abs. error for single-step Trap. : ", era(intval)
write(*,*) "Rel. error for single-step Trap. : ", err(intval)
write(*,*) "Integral from Trap. :", intval

write(11,*) log((b-a)), log(era(intval))

vals(1) = f(a)
vals(2) = f((a+b)/2.0)
vals(3) = f(b)

call simpint(1.75,3,valt,intval)

write(*,*) "Abs. error for single-step Simp. :", era(intval)
write(*,*) "Rel. error for single-step Simp, :", err(intval)
write(*,*) "Integral from Simp. :", intval

write(12,*) log((b-a)/2.0), log(era(intval))

! Deallocation of arrays
deallocate(vals)
deallocate(valt)

! Initial number of segments.
arrsize = 20

! Loop iterates the number of segments
do i = 1, 4, 1
	
	! Allocation and setting of dx
	allocate(valt(arrsize+1))
	dx = (b-a)/float(arrsize)

	!Finds values of f at xn
	do j = 1, arrsize+1, 1
		xi = a + float(j-1)*dx
		
		valt(j) = f(xi)
	
	end do
	
	! Calls integration subroutine
	call trapint(dx,arrsize+1,valt,intval)

	! Output
	write(*,51) arrsize, era(intval)
	write(*,52) arrsize, err(intval)
	write(*,*) "Integral from Trap. :", intval
	
	write(11,*) log(dx), log(era(intval))

	! Allocation and setting of dx
	allocate(vals(arrsize+1))
	dx = (b-a)/(float(arrsize))

	! Finds values of xn
	do j = 1, arrsize+1, 1
		xi = a + float(j-1)*dx
		
		vals(j) = f(xi)
	end do
	
	! Calls simpsons rule integration routine
	call simpint(dx,arrsize+1,vals,intval)
	
	! Output
	write(*,53) arrsize, era(intval)
	write(*,54) arrsize, err(intval)
	write(*,*) "Integral from Simp. :", intval
	
	write(12,*) log(dx), log(era(intval))
	
	! Iteration of number of segments
	arrsize = arrsize*2
	
	! Freeing of arrays
	deallocate(vals)
	deallocate(valt)
	
end do

end program

subroutine trapint(h,sz,arrin,intout)
! Subroutine containing algorithm for trapezoidal rule.

implicit none
	integer				:: sz
	real				:: h
	real, intent(in)	:: arrin(sz)
	real, intent(out)	:: intout

	integer				:: i
	real				:: partint
	
partint = 0.0
	
do i = 2, sz, 1

	partint = partint + h*(arrin(i) + arrin(i-1))/2.0

end do

	intout = partint
	
end subroutine

subroutine simpint(h,sz,arrin,intout)
! Subroutine containing algorithm for Simpsons rule.

implicit none
	integer				:: sz
	real				:: h
	real, intent(in)	:: arrin(sz)
	real, intent(out) 	:: intout
	
	integer				:: i
	real				:: partint
	
partint = 0.0

do i = 2, sz-2, 2

	partint = partint + 2.0*arrin(i) + 4.0*arrin(i+1)
	
end do

intout = (partint + arrin(1) + arrin(sz))*(h/3.0)

end subroutine
	
	
	
	
	
	
	
	


















