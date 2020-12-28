! Program to solve a 2D wave equation, with iteration in time.
! Spatial method will use a 5-point stencil, linearized, translated to a
! matrix-vector problem

! The iteration through time utilises the central time method.
! Both spatial and time methods together (CTCS) result in an explicit
! method of solution.

! Written for COMP 670, Homework 3

! Jon Parsons
! 10-11-19

program space

implicit none
  ! fn values, soln vector, coefficient matrix, inverse of coefficient matrix
  ! U vector - dim1 = space, dim2 = time
  ! A matrix - linearized stencil
  real*8, allocatable :: u_kwn(:,:), A_mat(:,:)
  ! Coefficient values
  real*8              :: c, R, err(3)
  integer             :: grd_max
  ! Error flag
  integer             :: all_err
  ! Discretization variables
  real*8              :: dx, dy, xi, yi, xf, yf
  real*8              :: dt, ti, tf
  ! Gridpoints in space and time
  integer             :: grdx, grdy, wrk_grd, wrk_grd2, grdt
  ! Looping integers
  integer             :: i

! maxima and minima of space
xi = 0.0
xf = 1.0
yi = 0.0
yf = 1.0
ti = 0.0
tf = 3.0
! Velocity of wave
c = 0.7
! Maximum number of gridpoints
grd_max = 160

! Smallest stepsize in space (dx = dy)
dx = (xf - xi)/real(grd_max-1,8)

! Set stepsize in time
dt = dx/(c*sqrt(2.0))

! Set number of steps in time
grdt = int((tf - ti)/dt)
grdt = grdt + 1 ! Plus one for time border

! Initial grid size
grdx = 20
grdy = 20

! Initialize error vector
err = 0.0

! Loop iterates the gridsize
m_iter: do i = 1, 3, 1

  ! Size of grid discounting the boundaries
  wrk_grd = (grdx - 2)
  wrk_grd2 = (wrk_grd)**2

  ! Determine grid size
  dx = (xf-xi)/float(grdx-1)
  dy = (yf-yi)/float(grdy-1)

  R = (c*dt/dx)**2 ! dx = dy

  ! Allocate arrays and intitalize
  allocate(A_mat(wrk_grd2,wrk_grd2), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 3. Exiting"
  A_mat = 0.0
  allocate(u_kwn(wrk_grd2,grdt+1), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 4. Exiting"
  u_kwn = 0.0
  write(*,*) "allocated"

  ! Fill A and initial values
  call ini_fill(wrk_grd,wrk_grd,grdt+1,A_mat,u_kwn,xi,yi,xf,yf,dx,dy,R)

  ! Solve the system
  call CTCS_solver(wrk_grd2,grdt+1,A_mat,u_kwn)

  ! Output data
  call outputs(wrk_grd,grdt+1,3,u_kwn,dx,i,err)

  ! Reset arrays
  deallocate(A_mat)
  deallocate(u_kwn)

  ! Update grid size
  grdx = grdx*2
  grdy = grdy*2

end do m_iter

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outputs(dim,dimt,err_dim,num,h,it_cn,err)
! Subroutine that outputs the results of the calculations

implicit none
  integer,intent(in)    :: dim, dimt, it_cn, err_dim
  real*8,intent(in)     :: num((dim)*(dim),dimt)
  real*8,intent(in)     :: h
  real*8,intent(inout)  :: err(err_dim)

  integer               :: i, j, k
  real*8                :: p, x_it, y_it

  character*20          :: file_num, err_file, file_time

51 format("num_grid",i0,".dat")
52 format(i0,"times.dat")

write(file_num,51) dim+2
err_file = "err_pts.dat"

! Open files
open(unit=15,file=trim(file_num),status="replace",position="append")

! Print values for plotting as a surface
do i = dim, 1, -1

  x_it = h*float(i)

  do j = dim, 1, -1

    y_it = h*float(j)

    write(15,*) x_it, y_it, num(j+((i-1)*dim),dimt)

    if ((i .eq. dim/2).and.(j .eq. dim/2)) then
      err(it_cn) = num(j+((i-1)*dim),dimt)
    end if

    end do
end do

close(15)

! Outputs the wave at all times. Uses 'X' as delimiter between times
if (dim+2 .eq. 20) then

  write(file_time,52) dim+2
  open(unit=16,file=file_time,status="replace",position="append")

  do k = 2, dimt, 1
    do i = dim, 1, -1

      x_it = h*float(i)

      do j = dim, 1, -1

        y_it = h*float(j)

        write(16,*) x_it, y_it, num(j+((i-1)*dim),k)

        end do
    end do


    write(16,*) "X"
  end do

end if

! Write point values
open(unit=17,file=trim(err_file),status="unknown",position="append")
write(17,*) err(it_cn)
close(17)


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CTCS_solver(dim,dimt,A,sln_mat)
! A subroutine which solves the Central-Time-Central-Space method. The
! takes the stencil matrix (A) and calls a multiplication routine against slices
! of the solution matrix, which holds the solutions at all times of interest.

implicit none
  integer,intent(in)    :: dim, dimt
  real*8,intent(in)     :: A(dim,dim)
  real*8,intent(inout)  :: sln_mat(dim,dimt)

  integer               :: i

! The time step being solved for is the i-th indice of the sln matrix
! making i-1 the current timestep, and i-2 the previous
do i = 3, dimt, 1

  call mv_mult(dim,A,sln_mat(:,i-1),sln_mat(:,i))
  sln_mat(:,i) = sln_mat(:,i) - sln_mat(:,i-2)

end do



end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ini_fill(dim1,dim2,grd_t,A,num,xi,yi,xf,yf,x_step,y_step,scale)
! Subroutine to fill the matrix and initial conditions of the vector.

implicit none
  ! Size of matrix/vector
  integer,intent(in)    :: dim1, dim2, grd_t
  ! Define matrix, vector within subroutine
  real*8,intent(inout)  :: A((dim1)*(dim2),(dim1)*(dim2))
  real*8,intent(inout)  :: num((dim1)*(dim2),grd_t)
  ! R, initial space, h
  real*8,intent(in)     :: xi, yi, xf, yf
  real*8,intent(in)     :: x_step, y_step, scale

  ! Function variables
  real*8                :: f, x, y, t, b, l, r
  real*8,parameter      :: pi = acos(-1.0) ! it's pi

  ! Looping variable
  integer               :: i, j, c
  ! Current value of x
  real*8                :: x_it, y_it

! Function, initial values
f(x,y) = sin(pi*x)*sin(pi*y)

! Boundary Conditions
! top
t(x,y) = 0.0*x*y
! bottom
b(x,y) = 0.0*x*y
! left
l(x,y) = 0.0*x*y
! right
r(x,y) = 0.0*x*y

! Fill initial values
do i = dim1, 1, -1

  x_it = xi + x_step*float(i)

  do j = dim1, 1, -1

    ! Fill RHS vector
    y_it = yi + y_step*float(j)

    ! Array stored S.T. y iterates fastest
    num(j+((i-1)*dim1),1) = f(x_it,y_it)
    num(j+((i-1)*dim1),2) = num(j+((i-1)*dim1),1)
    ! Boundary Conditions - Boundary values at x = 0,1 and y = 0,1
    ! Bottom (x = xi; y = yi)
    if (i .eq. 1) then
      num(j+((i-1)*dim1),1) = num(j+((i-1)*dim1),1) + b(x_it,yi)
      num(j+((i-1)*dim1),2) = num(j+((i-1)*dim1),1)

    end if
    ! Top (x = xi; y = yf)
    if (i .eq. (dim1)) then
      num(j+((i-1)*dim1),1) = num(j+((i-1)*dim1),1) + t(x_it,yf)
      num(j+((i-1)*dim1),2) = num(j+((i-1)*dim1),1)
    end if
    ! Left (x = xf; y = yi)
    if (j .eq. 1) then
      num(j+((i-1)*dim1),1) = num(j+((i-1)*dim1),1) + l(xi,y_it)
      num(j+((i-1)*dim1),2) = num(j+((i-1)*dim1),1)
    end if
    ! Right (x = xf; y = yi)
    if (j .eq. (dim1)) then
      num(j+((i-1)*dim1),1) = num(j+(dim1*(i-1)),1) + l(xf,y_it)
      num(j+((i-1)*dim1),2) = num(j+((i-1)*dim1),1)
    end if

  end do

end do

c = 1
! Fill main matrix
vert_loop: do i = 1, (dim1)*(dim2), 1

  ! Main diagonal
  A(i,i) = -4.0*scale + 2
  ! Super and sub diagonals
  if (i .lt. dim1*dim2) then
    A(i+1,i) = 1.0*(scale)
    A(i,i+1) = 1.0*(scale)
  end if
  ! Extremal super and sub diagonals
  if (i .gt. dim1) then
    A(i,i-(dim1)) = 1.0*(scale)
  end if
  if (i .le. (dim1*dim2-(dim2))) then
    A(i,i+dim1) = 1.0*(scale)
  end if
  ! Boundary conditions
  if (((i .gt. 1).and.(c .eq. (dim1))).and.(i .lt. dim1*dim2)) then
    A(i,i+1) = 0.0
    A(i+1,i) = 0.0
    c = 0
  end if

  c = c + 1

end do vert_loop

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
