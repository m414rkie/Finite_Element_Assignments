! Program to solve a 2D wave equation, no iteration in time.
! Method will use a 5-point stencil, linearized, translated to a
! matrix-vector problem

! Revised version.

! Previous iteration utlized a naive linearization which
! introduced errors into the system that were independent of step-size.

! This version implements a new method of tracking the (x,y) coordinates.

! Written for COMP 670, Homework 2

! Jon Parsons
! 9-24-19

program space

implicit none
  ! Holds the wave at all points (x,y). num => numerical  ana => analytical
  real*8, allocatable   ::  u_ana(:)
  ! fn values, soln vector, coefficient matrix, inverse of coefficient matrix
  real*8, allocatable   :: u_kwn(:), u_ukwn(:), A_mat(:,:), A_inv(:,:)
  ! Error flag
  integer             :: all_err
  ! Discretization variables
  real*8                :: dx, dy, xi, yi, xf, yf
  ! Gridpoints in space and time
  integer             :: grdx, grdy, wrk_grd, wrk_grd2
  ! Looping integers
  integer             :: i

! maxima and minima of space
xi = 0.0
xf = 1.0
yi = 0.0
yf = 1.0

! Initial grid size
grdx = 15
grdy = 15

! Loop iterates the gridsize
m_iter: do i = 1, 3, 1

  ! Size of grid discounting the boundaries
  wrk_grd = (grdx - 2)
  wrk_grd2 = (wrk_grd)**2

  ! Determine grid size
  dx = (xf-xi)/float(grdx-1)
  dy = (yf-yi)/float(grdy-1)

  ! Allocate arrays and intitalize
  allocate(u_ana(wrk_grd2), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 2. Exiting"
  u_ana = 0.0
  allocate(A_mat(wrk_grd2,wrk_grd2), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 3. Exiting"
  A_mat = 0.0
  allocate(A_inv(wrk_grd2,wrk_grd2), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 3. Exiting"
  A_inv = 0.0
  allocate(u_kwn(wrk_grd2), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 4. Exiting"
  u_kwn = 0.0
  allocate(u_ukwn(wrk_grd2), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 5. Exiting"
  u_ukwn = 0.0
  write(*,*) "allocated"

  ! Fill A and initial values
  call ini_fill(wrk_grd,wrk_grd,A_mat,u_ana,u_kwn,xi,yi,xf,yf,dx,dy)

  ! Invert stencil matrix
  call invert(A_mat,A_inv,wrk_grd2)
  ! Test inversion
  call mat_mul(wrk_grd2,A_mat,A_inv)
  ! Solve the system
  call mv_mult(wrk_grd2,A_inv,u_kwn,u_ukwn)

  ! Output data
  call outputs(wrk_grd,u_ana,u_ukwn,dx)

  ! Reset arrays
  deallocate(u_ana)
  deallocate(A_mat)
  deallocate(A_inv)
  deallocate(u_kwn)
  deallocate(u_ukwn)

  ! Update grid size
  grdx = grdx*2
  grdy = grdy*2

end do m_iter

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mat_mul(dim,A,B)
! Subroutine to multiply two square matrices for testing to see if inversion
! worked correctly

implicit none
  ! Dimension of matrices
  integer,intent(in)    :: dim
  ! Matrix 1, matrix 2
  real*8,intent(in)     :: A(dim,dim), B(dim,dim)

  ! Resulting matrix
  real*8                :: C(dim,dim)
  ! Holding variable, for use in finding difference
  real*8                :: val, sum
  ! Looping integers
  integer               :: i, j, k

! Initialize
C = 0.0

! Multiplication loops
do i = 1, dim, 1
  do j = 1, dim, 1
    val = 0.0
    do k = 1, dim, 1
      val = val + A(i,k)*B(k,j)
    end do
    C(i,j) = val
  end do
end do

! Sum resultant Matrix
sum = 0.0
do i = 1, dim, 1
  do j = 1, dim, 1
    sum = sum + c(i,j)
  end do
end do

! The resultant matrix should be the identity matrix with a summed value equal
! to the dimensions
write(*,*) "Error of Inversion:", sum - dim

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outputs(dim,ana,num,h)
! Subroutine that outputs the results of the calculations

implicit none
  integer,intent(in)    :: dim
  real*8,intent(in)     :: ana((dim)*(dim)), num((dim)*(dim)), h

  integer               :: i, j
  real*8                :: l2_norm, x_it, y_it

  character*20          :: file_num, file_kwn, err_file

51 format("num_grid",i0,".dat")
52 format("ana_grid",i0,".dat")
53 format("err_pts.dat")

write(file_num,51) dim+2
write(file_kwn,52) dim+2
write(err_file,53)

! Open files
open(unit=15,file=trim(file_num),status="replace",position="append")
open(unit=16,file=trim(file_kwn),status="replace",position="append")
open(unit=17,file=trim(err_file),status="unknown",position="append")

! Find the error
l2_norm = 0.0
do i = 1, (dim)*(dim), 1
  l2_norm = l2_norm + ((ana(i)-num(i))**2)
end do

l2_norm = h*sqrt(l2_norm)

write(17,*) h, l2_norm

! Print values for plotting as a surface
do i = dim, 1, -1

  x_it = h*float(i)

  do j = dim, 1, -1

    y_it = h*float(j)

    write(15,*) x_it, y_it, num(j+((i-1)*dim))
    write(16,*) x_it, y_it, ana(j+((i-1)*dim))

  end do
end do

close(15)
close(16)
close(17)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine invert(a,ainv,dime)
! Inversion subroutine. Finds the inverse of matrix A column-wise by
! solving [A:I] -> [I:A^-1]

implicit none
  ! Dimension of matrix
	integer				       :: dime
  ! Input and resulting inverse
	real*8				      :: a(dime,dime), ainv(dime,dime)
  ! Looping integers and size of working matrix
	integer				      :: i,j, dub
  ! Working matrix
	real*8,allocatable	:: work(:,:)

dub = 2*dime

allocate(work(dime,dub))

! Fill working matrix
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
	! make appropriate division to row in RHS
	work(i,i+1:dub) = work(i,i+1:dub)/work(i,i)
	! Set main diagonal
	work(i,i) = 1.0
	do j = 1, dime, 1
		if (j .eq. i) then
			cycle
		end if
		! Remainder of row work
		work(j,i+1:dub) = work(j,i+1:dub)-work(i,i+1:dub)*work(j,i)
		work(j,i) = 0.0
	end do
end do

! Fill outgoing matrix by column
do i = 1, dime, 1
	ainv(:,i) = work(:,i+dime)
end do

deallocate(work)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ini_fill(dim1,dim2,A,ana,num,xi,yi,xf,yf,x_step,y_step)
! Subroutine to fill the matrix and initial conditions of the vector.

implicit none
  ! Size of matrix/vector
  integer,intent(in)    :: dim1, dim2
  ! Define matrix, vector within subroutine
  real*8,intent(inout)  :: A((dim1)*(dim2),(dim1)*(dim2))
  real*8,intent(inout)  :: num((dim1)*(dim2)), ana(dim1*dim2)
  ! R, initial space, h
  real*8,intent(in)     :: xi, yi, x_step, y_step, xf, yf

  ! Function variables
  real*8                :: g, f, x, y, t, b, l, r

  ! Looping variable
  integer               :: i, j, c
  ! Current value of x
  real*8                :: x_it, y_it
  ! Inverse of both 'steps'
  real*8                :: h2, h2_inv

! Function, analytical solution
!g(x,y) = x+y
g(x,y) = x*y*(x-1.0)*(y-1.0)*exp(x-y)

! Function, RHS values
f(x,y) = -2.0*x*(y-1.0)*(y-2.0*x+x*y+2.0)*exp(x-y)

! Boundary Conditions
! top
t(x,y) = 0.0*x*y
! bottom
b(x,y) = 0.0*x*y
! left
l(x,y) = 0.0*x*y
! right
r(x,y) = 0.0*x*y

h2 = (x_step*y_step)
h2_inv = 1.0/h2

! Fill RHS
do i = dim1, 1, -1

  x_it = xi + x_step*float(i)

  do j = dim1, 1, -1

    ! Fill RHS vector
    y_it = yi + y_step*float(j)

    ! Array stored S.T. y iterates fastest
    num(j+((i-1)*dim1)) = f(x_it,y_it)
    ! Boundary Conditions - Boundary values at x = 0,1 and y = 0,1
    ! Bottom (x = xi; y = 0)
    if (i .eq. 1) then
      num(j+((i-1)*dim1)) = num(j+((i-1)*dim1)) + b(x_it,yi)
    end if
    ! Top (x = xi; y = 1)
    if (i .eq. (dim1)) then
      num(j+((i-1)*dim1)) = num(j+((i-1)*dim1)) + t(x_it,yf)
    end if
    ! Left (x = 0; y = yi)
    if (j .eq. 1) then
      num(j+((i-1)*dim1)) = num(j+((i-1)*dim1)) + l(xi,y_it)
    end if
    ! Right (x = 1; y = yi)
    if (j .eq. (dim1)) then
      num(j+((i-1)*dim1)) = num(j+(dim1*(i-1))) + l(xf,y_it)
    end if

  end do

end do

! Fill analytic vector
do i = dim1, 1, -1

  x_it = xi + x_step*float(i)

  do j = dim1, 1, -1

    y_it = yi + y_step*float(j)

    ana(j+((i-1)*dim1)) = g(x_it,y_it)

  end do

end do
c = 1
! Fill main matrix
vert_loop: do i = 1, (dim1)*(dim2), 1

  ! Main diagonal
  A(i,i) = -4.0
  ! Super and sub diagonals
  if (i .lt. dim1*dim2) then
    A(i+1,i) = 1.0
    A(i,i+1) = 1.0
  end if
  ! Extremal super and sub diagonals
  if (i .gt. dim1) then
    A(i,i-(dim1)) = 1.0
  end if
  if (i .le. (dim1*dim2-(dim2))) then
    A(i,i+dim1) = 1.0
  end if
  ! Boundary conditions
  if (((i .gt. 1).and.(c .eq. (dim1))).and.(i .lt. dim1*dim2)) then
    A(i,i+1) = 0.0
    A(i+1,i) = 0.0
    c = 0
  end if

  c = c + 1

end do vert_loop

! Adjust to stepsize
A = -h2_inv*A

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mv_mult(dim,A,V,O)
! Subroutine to multiply a square matrix and vector of corresponding size.

implicit none
  ! Dimension of matrix/vector
  integer,intent(in)   :: dim
  ! Define matrix and vector
  real*8,intent(in)    :: A(dim,dim), V(dim)
  ! Define output vector
  real*8,intent(out)   :: O(dim)

  ! Looping integers
  integer              :: i, j

! Initialize
O = 0.0

! Row iteration
do i = 1, dim, 1
  ! Column iteration
  do j = 1, dim, 1
      ! The calculation
      O(i) = O(i) + A(i,j)*V(j)

  end do

end do

end subroutine
