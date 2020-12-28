program difference

!	Written for use in completing homework two for COMP521. 
! Calculates the derivative of two functions using the three
! first order finite differences methods as well as a second 
! derivative finite difference method. Program iterates dx
! with outer loop. Files are renamed on the fly for output
! into a directory in the same location as the executable.
! Derivative methods are contained in subroutines.
!
! Author: Jon Parsons, 9-29-2018

implicit none
! Common Variables
integer				:: i, j, k							! Looping integers
real				:: f, g, x							! Function variables
real				:: exactff, exactfg, exactsff		
real				:: exactf, exactg, exactsf		
real, allocatable	:: arr1(:), arr2(:)					! Array variables
real, allocatable	:: dfor1(:,:), dback1(:,:)
real, allocatable	:: dcent1(:,:)
real, allocatable	:: dfor2(:,:), dback2(:,:)
real, allocatable	:: dcent2(:,:)
character(50)		:: file1, file2, file3				! Filenames
character(50)		:: filepath							! Directory
! Differentiation Variables
real				:: dx								! Stepsize variable
real				:: dom, a, b						! Domain size and endpoints
integer				:: numstep, un						! Count variables
real				:: err, rang
! Calculation Variables
real				:: xi								! Current x value

! Function Statements
f(x) = sin(2.0*x)
g(x) = log(3.0*x + 1.0)
exactff(x) = 2.0*cos(2.0*x)
exactfg(x) = 3.0/(3.0*x + 1.0)
exactsff(x) = -4.0*sin(2.0*x) 

! Format Strings
50 format ("Output/",f6.4,"foreqn1.dat")
51 format ("Output/",f6.4,"foreqn2.dat")
52 format ("Output/",f6.4,"backeqn1.dat")
53 format ("Output/",f6.4,"backeqn2.dat")
54 format ("Output/",f6.4,"centeqn1.dat")
55 format ("Output/",f6.4,"centeqn2.dat")

! Checks for filepath
filepath = "~/Desktop/Comp521/HW2_Differentiation/Output"
call dircheck(filepath)

! Initializations
a = 0.5
b = 2.5
dom = (b - a)
dx = 0.1

! Outer loop, iterates dx
do k = 1, 4, 1

	numstep = int(dom/dx)
	
	! Allocation statements
	allocate(arr1(numstep))    ; allocate(arr2(numstep))
	allocate(dfor1(3,numstep)) ; allocate(dback1(3,numstep))
	allocate(dcent1(3,numstep))
	allocate(dfor2(3,numstep)) ; allocate(dback2(3,numstep))
	allocate(dcent2(3,numstep))

	
	! Calculation statements
	! Evaluates functions
	do i = 1, numstep, 1

		xi = dx*float(i) + a

		dback1(1,i) = xi ; dback1(3,i) = exactff(xi)
		dfor1(1,i) = xi  ; dfor1(3,i) = exactff(xi)
		dcent1(1,i) = xi ; dcent1(3,i) = exactff(xi)
		dfor2(1,i) = xi	 ; dfor2(3,i) = exactfg(xi)
		dback2(1,i) = xi ; dback2(3,i) = exactfg(xi)
		dcent2(1,i) = xi ; dcent2(3,i) = exactfg(xi)
		arr1(i) = f(xi)
		arr2(i) = g(xi)
	
	end do
		
	! Calls diferentiation subroutines
	call backward(arr1,numstep,dx,dback1)
	call backward(arr2,numstep,dx,dback2)
	call forward(arr1,numstep,dx,dfor1)
	call forward(arr2,numstep,dx,dfor2)
	call central(arr1,numstep,dx,dcent1)
	call central(arr2,numstep,dx,dcent2)

	write(file1,50) dx
	write(file2,51) dx

	open(unit=21,file=trim(file1),status="replace",position="append")
	open(unit=22,file=trim(file2),status="replace",position="append")
			
	write(21,*) "x	finite	exact"
	write(22,*) "x	finite	exact"
	
	! Write Loops
	do j = 1, numstep-1, 1
		write(21,*) dfor1(1,j), dfor1(2,j), dfor1(3,j)
		write(22,*) dfor2(1,j), dfor2(2,j), dfor2(3,j)
	end do
	
	close(21)
	close(22)
	
	
	write(file1,52) dx
	write(file2,53) dx

	open(unit=21,file=trim(file1),status="replace",position="append")
	open(unit=22,file=trim(file2),status="replace",position="append")
	
	write(21,*) "x	finite	exact"
	write(22,*) "x	finite	exact"
	
	do j = 2, numstep, 1
		write(21,*) dback1(1,j), dback1(2,j), dback1(3,j)
		write(22,*) dback2(1,j), dback2(2,j), dback2(3,j)
	end do
	
	close(21)
	close(22)

	
	write(file1,54) dx
	write(file2,55) dx

	open(unit=21,file=trim(file1),status="replace",position="append")
	open(unit=22,file=trim(file2),status="replace",position="append")
	
	write(21,*) "x	finite	exact"
	write(22,*) "x	finite	exact"
	
	do j = 2, numstep-1, 1
		write(21,*) dcent1(1,j), dcent1(2,j), dcent1(3,j)
		write(22,*) dcent2(1,j), dcent2(2,j), dcent2(3,j)
	end do
	
	close(21)
	close(22)
	
	! Updates dx
	dx = dx*0.5

	deallocate(arr1)  ; deallocate(arr2)
	deallocate(dfor1) ; deallocate(dback1)
	deallocate(dcent1)
	deallocate(dfor2) ; deallocate(dback2)
	deallocate(dcent2)

end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here begins Code for Problem 2
! Similar structure are part one

56 format ("Output/",f6.4,"secder.dat")
open(unit=28,file="Output/Errplot.dat",status="replace",position="append")

dx = 0.1
un = 1
rang = 0.01

do k = 1, 4, 1

	numstep = int(dom/dx)	
	write(file3,56) dx
	
	allocate(arr1(numstep))	
	allocate(dfor1(3,numstep))

	do i = 1, numstep, 1

	
	xi = dx*float(i) + a
	arr1(i) = f(xi)
	dfor1(1,i) = xi ; dfor1(3,i) = exactsff(xi)

	
	if (abs(xi-1.0) .lt. rang) then
		un = i
	end if

	end do
	
	call secder(arr1,numstep,dx,dfor1)
	
	open(unit=22,file=trim(file3),status="replace",position="append")
	write(22,*) "x	finite	exact"
	
	! Error and writing loop.
	! Error values taken from x=1
	err = 0.0

	do j = 3, numstep-2, 1
		write(22,*) dfor1(1,j), dfor1(2,j), dfor1(3,j)
		err = dfor1(un,2) - dfor1(un,3)	
	end do


	close(22)
	
	write(28,*) log(dx), log(abs(err))
	
	dx = dx*0.5
	
	deallocate(arr1)	
	deallocate(dfor1)
	
end do

close(28)

end program

subroutine forward(arrin,sz,step,out)
! Forward differentiation subroutine
implicit none
	integer,intent(in)		:: sz
	real, intent(in)		:: arrin(sz)
	real,intent(in)			:: step
	real,intent(out)		:: out(3,sz)
	
	integer			:: i

do i = 1, sz-1, 1
	out(2,i) = (arrin(i+1) - arrin(i))/step
end do

end subroutine

subroutine backward(arrin,sz,step,out)
! Backwards differentiation subroutine
implicit none
	integer,intent(in)		:: sz
	real, intent(in)		:: arrin(sz)
	real,intent(in)			:: step
	real,intent(out)		:: out(3,sz)
	
	integer			:: i

do i = 2, sz, 1
	out(2,i) = (arrin(i) - arrin(i-1))/step
end do

end subroutine

subroutine central(arrin,sz,step,out)
! Central differentiation subroutine
implicit none
	integer,intent(in)		:: sz
	real, intent(in)		:: arrin(sz)
	real,intent(in)			:: step
	real,intent(out)		:: out(3,sz)
	
	integer			:: i

do i = 2, sz-1, 1
	out(2,i) = (arrin(i+1) - arrin(i-1))/(2.0*step)
end do

end subroutine

subroutine secder(arrin,sz,step,out)
! Second derivitive subroutine
implicit none
	integer,intent(in)		:: sz
	real,intent(in)			:: arrin(sz)
	real,intent(in)			:: step
	real,intent(out)		:: out(3,sz)
	
	integer					:: i

do i = 3, sz-2, 1	
	out(2,i) = (1.0/(12.0*(step**2)))*(-arrin(i-2)+16.0*arrin(i-1)-30.0*arrin(i)+16.0*arrin(i+1)-arrin(i+2))
end do

end subroutine
	
subroutine dircheck(path)

implicit none
	character(len=*)		:: path
	character(len=100)		:: makepath
	logical					:: check
	
	
	inquire(file=trim(path)//'/.', exist=check)
	
	if (check .eqv. .false.) then
		makepath = "mkdir -p "//trim(path)
		call system(makepath)
	end if

end subroutine
	