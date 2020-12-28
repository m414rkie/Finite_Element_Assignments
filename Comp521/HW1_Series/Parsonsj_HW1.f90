program series

! Written for use in completing homework one for COMP521
! Author: Jon Parsons, 9-13-2018

implicit none
! Common variables
integer			:: i
real			:: x, y
! Variables for question 1.A
real			:: sa, apartial
real, parameter :: pi = acos(-1.0)
! Variables for question 1.B
real			:: sb, bpartial
! Variables for question 1.C
real			:: sc, cpartial
! Variables for question 2
real			:: exact, t2, t5
real			:: norm2, norm5
real			:: norm2quad(1:500), norm5quad(1:500)

! Declaring functions

! 	Defining functions for question 1
! Series A given as Stirling's approximation. Factorials above approx. 20 require
! unreasonable amounts of memory. Approximation: ln(n!) = nln(n)-n+ln(sqrt(2pin))

sa(x) = (x*log(x)-x+log(sqrt(2.0*pi*x)))**2/(2.0*x*log(2.0*x)-2.0*x+log(sqrt(4.0*pi*x)))

sb(x) = (3**(2.0*x))/(10**x)

sc(x) = exp(-(x**2))

!	Defining functions for question 2, Taylor series centered on 0.

exact(x) = log(3.0*x + 1.0)

! First three terms of Taylor Series
t2(x) = log(1.0) + (3.0*x/(1.0))

! Terms four and five
t5(x) = ((-9.0*(x**2))/(2.0)*((1.0)**2)) + (54.0*(x**3))/((6.0)*(1.0)**3) + (-486.0*(x**4))/((24.0)*((1.0)**4))

! Opening files

open(unit=11,file="HW1_series.dat",status="replace",position="append")
write(11,*) "n	A		B		C"

open(unit=12,file="HW1_taylor.dat",status="replace",position="append")
write(12,*) "n	Exact	t2		t5"

! Initializing variables
apartial = 0.0
bpartial = 0.0
cpartial = 0.0
norm2quad = 0.0
norm5quad = 0.0

! Operations for question 1

do i = 1, 50, 1

	apartial = apartial + sa(real(i))
	bpartial = bpartial + sb(real(i-1))
	cpartial = cpartial + sc(real(i))
		
	write(11,*) i, apartial, bpartial, cpartial

end do


! Operations for question 2

do i = 1, 500, 1
	
	y = float(i)*0.001
	
	write(12,*) y, exact(y), t2(y), t2(y)+t5(y)
	norm2quad(i) = (exact(y) - t2(y))**2
	norm5quad(i) = (exact(y) - (t2(y)+t5(y)))**2
	
end do

call boolequad(norm2quad,norm2,500)
write(*,*) norm2

call boolequad(norm5quad,norm5,500)
write(*,*) norm5

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine boolequad(arr,booleint,upto)

! Integration subroutine, uses boolean quadrature method (5-pt)

implicit none
	integer					  :: upto		! Array input size
	real, dimension(upto)	  :: arr		! Dummy name for array
	real, intent(out)		  :: booleint	! Final value for this subroutine
	integer					  :: i			! Looping integer

! Initialization
booleint = 0.0

! Loop for composite Boole's rule. Weights included.
do i = 1, upto, 5
	
	booleint = booleint + 0.001*(2.0/45.0)*(7.0*arr(i+4)+32.0*arr(i+3)+12.0*arr(i+2)+32.0*arr(i+1)+7.0*arr(i))

end do

end subroutine

