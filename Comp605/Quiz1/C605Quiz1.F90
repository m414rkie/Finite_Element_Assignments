! Code written for quiz one, COMP 605
! Code will take two matrices and multiply them
! using collective communications in MPI

! Method calculates rows independently. Method
! of parallelization cannot handle odd numbers
! of cores at this point. 

! AUthor: Jon Parsons
! 3-12-2019

program matmulcoll

use mpi

implicit none
	integer, parameter	:: n = 8 ! size of matrices, square
	real				:: a(n,n), b(n,n) ! Matrices to be multiplied
	integer				:: h ! Stepsize for passing rows/columns to subroutines
	integer				:: i_i, i_f ! Initial and final values of loops for multiplication
	real				:: c(n,n), loc_c(n,n) ! Global and local results multiplication
	
	integer				:: world_size	! Holds number of processors available
	integer				:: mpi_err, ierr ! MPI error variable
	integer				:: local_num ! Holds local address of processor
	
	integer				:: i, j ! looping integers
	

! Fill matrices a, b, and c
call matfill(n,a,b)
c = 0.0


! Initialize MPI
call MPI_INIT (mpi_err)
call MPI_COMM_SIZE (MPI_COMM_WORLD,world_size, mpi_err)
call MPI_COMM_RANK (MPI_COMM_WORLD,local_num, mpi_err)

! Broadcast global information
call MPI_BCAST(a,n,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(b,n,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

! determine number of rows per process
h = n/world_size

! Initialize
loc_c = 0.0

! Determine local endpoints
i_i = 1 + local_num*h
i_f = i_i + h - 1

call matmult(n,loc_c,a,b,i_i,i_f,1,n)

! Collect results
call MPI_REDUCE(loc_c,c,n*n,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,mpi_err)

! Output to user
if (local_num .eq. 0) then
	write(*,*) "Program Complete."
	write(*,*) "Matrix A:"
	do i = 1, n, 1
		do j = 1, n, 1
			write(*,'(1f4.1, " ")',advance="no") a(i,j)
		end do
		write(*,*)
	end do
	
	write(*,*) 
	write(*,*) "Matrix B:"
	do i = 1, n, 1
		do j = 1, n, 1
			write(*,'(1f4.1, " ")',advance="no") b(i,j)
		end do
		write(*,*)
	end do
	write(*,*) 
	write(*,*) "Results of multiplication:"
	do i = 1, n, 1
		do j = 1, n, 1
			write(*,'(1f6.1, " ")',advance="no") c(i,j)
		end do
		write(*,*)
	end do
end if

call MPI_FINALIZE(mpi_err)
	
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matfill(dime,matA,matB)
! Subroutine fills initial matrix values 

implicit none
	integer		:: dime
	real		:: matA(dime,dime), matB(dime,dime)
	
	integer		:: i, j

matA = 0.0
matB = 0.0

do i = 1, dime, 1
	do j = 1, dime, 1
		matA(i,j) = 2.0*float(i) + float(j)
		matB(i,j) = 3.0*float(j) - float(i)
	end do
end do

end subroutine
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
subroutine matmult(dime,matC,matA,matB,i1,i2,j1,j2)
! subroutine multiplies the matrices together

implicit none
	integer		:: dime ! dimensions
	real		:: matC(dime,dime), matA(dime,dime), matB(dime,dime) ! matrices
	integer		:: i1,i2,j1,j2 ! initial and final endpoints of loops
	
	
	integer		:: i, j, k ! looping integers
		
matC = 0.0

do i = i1, i2, 1
	do j = j1, j2, 1
		do k = 1, dime, 1
			matC(i,j) = matC(i,j) + matA(i,k)*matB(k,j)
		end do
	end do
end do
	
	
end subroutine
		
	
		
		
		
		
		
		




