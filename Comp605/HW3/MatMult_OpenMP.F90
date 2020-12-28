! Code for HW5, COMP605.

! Program will multiply two matrices together using an
! OpenMP parallelization scheme.

! Compile as : gfortran -o MatMult.x MatMult_OpenMP.F90 -fopenmp

! Be sure to set number of threads by command line
! EX: export OMP_NUM_THREADS=#

! Author: Jon Parsons

program MatMul_openMP

use omp_lib

implicit none
	integer,parameter	:: N=16
	real				:: A(N,N), B(N,N) ! The matrices to be multiplied together
	real				:: C(N,N)	! The resultant matrix
	integer				:: i, j, k	! Looping integers
	integer				:: NumThreads ! Number of threads to be used, user input
	integer,parameter	:: dp = kind(1.d0)
	real (kind=dp)		:: walltime ! Holds the timing
	
write(*,*) "Please enter the number of threads to use:"
read(*,*) NumThreads

! Fill initial matrices (Serial)
do i = 1, N, 1
	do j = 1, N, 1
		A(i,j) = float(i+j)
		B(i,j) = float(i-j)
	end do
end do

! Initialize C (Serial)
C = 0.0

! Initializa OpenMP variable

!$omp parallel shared (A,B,C) private (i,j,k)

! Define number of threads
call omp_set_num_threads(NumThreads)

! Timing calls
walltime = omp_get_wtime()

! Multiplication loop
!$omp do
do i = 1, N, 1
	do j = 1, N, 1
		do k = 1, N, 1
			C(i,j) = c(i,j) + A(i,k)*B(k,j)
		end do
	end do
end do
!$omp end do

!$omp end parallel

! Finaliza timing
walltime = omp_get_wtime() - walltime
write(*,*) "Elapsed time (ms): ", walltime*1000

! Display C
!do i = 1, N, 1
!	do j = 1, N, 1
!		write(*,'(f10.2)',advance="NO") C(i,j)
!	end do
!	write(*,*)
!end do

end program


