program CG

! Program written for problem one of the second homework assignment in Comp605.
! Program solves a system of equations using the conjugate gradient method.
! The user will enter the size of a square matrix A. All vectors will be
! sized appropriately to ensure compatibility with the matrix.

! Author: Jon Parsons
! Date: 2-27-2019

implicit none
	integer				:: i, j, k ! Looping integers and k for iteration counting.
	integer				:: DimA 	! Array size
	real				:: alpha, beta ! Correction variables
	real				:: diff, err	! Loop exit condition variables
	real,allocatable	:: matA(:,:), vecB(:)	! Matrix A and vector B
	real,allocatable	:: matM(:)		! One dimensional to save memory. Jacobi preconditioning
										! only uses diagonal
	real,allocatable	:: vecR(:), vecR_p1(:) ! Residual vectors
	real,allocatable	:: vecP(:), vecP_p1(:) ! Correction vectors
	real,allocatable	:: vecX(:), vecX_p1(:) ! Variable Vectors
	real,allocatable	:: vecAp(:), vecAp_p1(:) ! Temporary holding vectors for use in finding 
													! alpha, beta
	real,allocatable	:: vecZ(:), vecZ_p1(:) ! For use with preconditioning
	
	integer				:: preflag	! Flag, if !=0 preconditioning is used
	real				:: rr, r1r1, pAp ! Hold values for finding alpha, beta.
	
	character*30		:: filename ! Filename for output

! User input
write(*,*) "Please enter the size of the square matrix. Each vector will be sized as appropriate."
read(*,*) DimA

write(*,*) "Preconditioning? (0 for no, 1 for yes)"
read(*,*) preflag

! Allocation statements
allocate(matA(DimA,DimA))
allocate(vecB(DimA))
allocate(vecR(DimA))
allocate(vecR_p1(DimA))
allocate(vecP(DimA))
allocate(vecP_p1(DimA))
allocate(vecX(DimA))
allocate(vecX_p1(DimA))
allocate(vecAp(DimA))
allocate(vecAp_p1(DimA))
allocate(matM(DimA))
allocate(vecZ(DimA))
allocate(VecZ_p1(DimA))

! Initialize the error allowed and iteration number
err = 0.000001
k = 0

! Filename logic, with preconditioning _pref is appended.
if (preflag .eq. 0) then
	write(filename,*) "convergence.dat"
else
	write(filename,*) "convergence_pref.dat"
end if

! Open file and write headers to data
open(unit=15,file=trim(filename),status="replace",position="append")
write(15,*) "Iteration	Error"

! Fill the arrays
matA = 0.0

do i = 1, DimA, 1
	
	matA(i,i) = 0.5 + sqrt(float(i))	
	matM(i) = 1.0/matA(i,i)
	
	do j = 1, DimA, 1
		
		if (abs(i-j) .eq. 100) then
			matA(i,j) = 1.0
		end if

	end do
	
end do

vecB = 1.0	
vecX = 0.0
vecX_p1 = 0.0

! Begin CG method
! This section before the loop checks if the current values of x are sufficient
! to exit
call matVec(DimA,matA,vecX,vecAp)

vecR = vecB - vecAp
if (preflag .eq. 0 ) then
	vecP = vecR
else
	do i = 1, DimA, 1
		vecZ(i) = matM(i)*vecR(i)
	end do
	vecP = vecZ
end if 

diff = abs(sum(vecR))

! Main working loop
do while (diff .gt. err)
		
		k = k + 1 ! Update iteration
		
		! Logic for including/disallowing preconditioning
		if (preflag .eq. 0) then
			call vecVec(DimA,vecR,vecR,rr)
		else
			call vecVec(DimA,vecR,vecZ,rr)
		end if
		
		call matVec(DimA,matA,vecP,vecAp)
		call vecVec(DimA,vecP,vecAp,pAp)
		
		alpha = rr/pAp

		vecX_p1 = vecX + alpha*vecP
		vecR_p1 = vecR - alpha*vecAp
		
		diff = abs(sum(vecR_p1))
		
		! Exit logic
		if (diff .lt. err) then
			write(15,*) k, log(diff)
			write(*,*) k, diff
			goto 101
		end if 
		
		! Logic for preconditioning
		if (preflag .eq. 1) then
			
			do i = 1, DimA, 1
				vecZ_p1(i) = matM(i)*vecR_p1(i)
			end do
			
			call vecVec(DimA,vecZ_p1,vecR_p1,r1r1)
			call vecVec(DimA,vecZ,vecR,rr)
			beta = r1r1/rr
			vecP_p1 = vecZ_p1 + beta*vecP
			
		else
		
			call vecVec(DimA,vecR_p1,vecR_p1,r1r1)
			beta = r1r1/rr
			vecP_p1 = vecR_p1 + beta*vecP
		
		end if 
		
		! Write current iteration and error to file.
		write(15,*) k, log(diff)
		! Update vectors for the next iteration
		vecX = vecX_p1
		vecP = vecP_p1
		vecR = vecR_p1
		vecZ = vecZ_p1
		
end do

101 write(*,*) "Solution Found"

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine vecVec(dime,vecA,vecB,partial)
! This subroutine multiplies two vectors together to find a scalar.

implicit none
	integer,intent(in)	:: dime
	real,intent(in)		:: vecA(dime), vecB(dime)
	real,intent(out)	:: partial
	
	integer				:: i
	
partial = 0.0	
	
do i = 1, dime, 1
	
	partial = partial + vecA(i)*vecB(i)

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matVec(dime,mat,vecA,vecB)
! This subroutine multiplies a square matrix with a vector to output a vector.

implicit none
	integer,intent(in)		:: dime
	real,intent(in)			:: mat(dime,dime)
	real,intent(in)			:: vecA(dime)
	real,intent(out)		:: vecB(dime)
	
	integer					:: i, j
	real					:: partial

vecB = 0.0
	
do i = 1, dime, 1
	
	do j = 1, dime, 1
		
		vecB(i) = vecB(i) + mat(i,j)*vecA(j)
	
	end do
	
end do

end subroutine








