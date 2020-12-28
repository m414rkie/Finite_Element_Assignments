! Program to solve the 1D heat equation using a Forward-Time-Centered-Space
! scheme. Scheme to be implemented using a matrix-vector method.

! Written for COMP 670, Homework 1

! Jon Parsons
! 9-10-19

program heat

implicit none
  ! Holds the wave at all points (x,t). num => numerical : ana => analytical
  real, allocatable   :: u_num(:,:), u_ana(:,:)
  real, allocatable   :: u_vec_1(:), u_vec_2(:), A_mat(:,:)
  ! Error flag
  integer             :: all_err
  ! Discretization variables
  real                :: dx, dt, xi, ti, xf, tf
  ! Gridpoints in space and time
  integer             :: grdx, grdt
  ! FTCS scheme variables
  real                :: r
  ! Analytical function variables
  real                :: g, x, y, x_a, t_a
  real,parameter      :: pi = acos(-1.0) ! pi
  ! Value of thermal diffusivity
  real                :: alpha = 0.25

  ! Looping integers
  integer             :: i, j, k

! Function, analytical solution
g(x,y) = sin(pi*x*0.5)*exp(-alpha*pi*pi*y*0.25)

! maxima and minima of time and space
ti = 0.0
tf = 2.0
xi = 0.0
xf = 2.0

! Initial grid size
grdx = 10

m_iter: do i = 1, 3, 1
  ! Determine grid size
  dx = (xf-xi)/float(grdx-1)
  dt = (dx*dx)/(2.0*abs(alpha))
  grdt = int((tf-ti)/dt)
  write(*,*) "Points:", grdx

  ! Allocate arrays and intitalize
  allocate(u_num(grdt,grdx), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 1. Exiting"
  u_num = 0.0
  allocate(u_ana(grdt,grdx), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 2. Exiting"
  u_ana = 0.0
  allocate(A_mat(grdx,grdx), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 3. Exiting"
  A_mat = 0.0
  allocate(u_vec_1(grdx), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 4. Exiting"
  u_vec_1 = 0.0
  allocate(u_vec_2(grdx), stat=all_err)
  if (all_err .ne. 0) stop "Failed to allocate arrays 5. Exiting"
  u_vec_2 = 0.0
  write(*,*) "allocated"

  ! Calculate r, ensure stability conditions
  r = alpha*dt/(dx*dx)

  ! Fill A and initial values
  call ini_fill(grdx,A_mat,u_vec_1,r,xi,dx)

  ! Begin iteration through time
  time_loop: do j = 1, grdt, 1
    ! Numerical sln portion
    call mv_mult(grdx,A_mat,u_vec_1,u_vec_2)
    ! Save and update values
    u_num(j,:) = u_vec_2
    u_vec_1 = u_vec_2
    ! Analytical sln portion
    do k = 1, grdx + 1, 1

      x_a = xi + float(k-1)*dx
      t_a = ti + float(j)*dt
      ! Logic to ensure array boundaries are not violated
      if (k .lt. (grdx+1)) then
        u_ana(j,k) = g(x_a,t_a)
      end if
    end do

  end do time_loop

  ! Output data
  call outputs(grdt,grdx,u_num,u_ana,dt,dx)

  ! Reset arrays
  deallocate(u_num)
  deallocate(u_ana)
  deallocate(A_mat)
  deallocate(u_vec_1)
  deallocate(u_vec_2)

  ! Update grid size
  grdx = grdx*2

end do m_iter

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine outputs(dim1,dim2,num,ana,t_step,x_step)
! Subroutine to analize the outputs using matrices containing the full spatial
! and time results

implicit none
  integer,intent(in)    :: dim1, dim2 ! time, space
  ! Holds numerical and analytical solutions
  real,intent(in)       :: num(dim1,dim2),ana(dim1,dim2)
  ! h and k respectively
  real,intent(in)       :: t_step, x_step
  ! Looping integers
  integer               :: i, j
  ! Hold errors
  real                  :: l2_norm, pt_err
  ! File names
  character*30          :: file_err, file_val, file_ana

! Character variable blanks
51 format("err_grid.dat")
52 format("num_grid",i2,".dat")
53 format("ana_grid",i2,".dat")
! Finalized the character variables
write(file_err,51)
write(file_val,52) dim2
write(file_ana,53) dim2

! Open files
open(unit=15,file=trim(file_err),status="unknown",position="append")
open(unit=16,file=trim(file_val),status="unknown",position="append")
open(unit=17,file=trim(file_ana),status="unknown",position="append")

! Initialize
l2_norm = 0.0
pt_err = 0.0

! Define error on at tf and at the middle of space
pt_err = (ana(dim1,dim2/2) - num(dim1,dim2/2))

time_loop: do i = 1, dim1, 1

  space_loop: do j = 1, dim2, 1

    ! Update Euclidean norm
    l2_norm = l2_norm + t_step*x_step*(ana(i,j)-num(i,j))**2

    ! write out both solutions for plotting
    write(16,*) (t_step*float(i)), (x_step*float(j-1)), num(i,j)
    write(17,*) (t_step*float(i)), (x_step*float(j-1)), ana(i,j)

  end do space_loop

end do time_loop

! Finalize value for l2 norm
l2_norm = sqrt(l2_norm)

! Write to file: gridsize, point error, ln(h), ln(L2)
write(15,*) dim2, abs(pt_err), log(x_step), log(l2_norm)

! Close files
close(15)
close(16)
close(17)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ini_fill(dim,A,V,val,xi,x_step)
! Subroutine to fill the matrix and initial conditions of the vector.

implicit none
  ! Size of matrix/vector
  integer,intent(in)    :: dim
  ! Define matrix, vector within subroutine
  real,intent(inout)    :: A(dim,dim),V(dim)
  ! R, initial space, h
  real,intent(in)       :: val, xi, x_step

  ! Function variables
  real                  :: f, x
  real,parameter        :: pi = acos(-1.0) ! pi

  ! Looping variable
  integer               :: i
  ! Current value of x
  real                  :: x_it

! Function, initial conditions
f(x) = sin(pi*x*0.5)

! Working loop, all assignments done within
do i = 1, dim+1, 1

  ! Fill initial spatial values
  x_it = xi + float(i-1)*x_step
  if (i .lt. (dim+1)) then
    V(i) = f(x_it)
  end if

  ! Fill main diagonal of matrix
  if (((i .eq. 1).or.(i .eq. dim-1).and.(i .le. (dim)))) then
    A(i,i) = 1.0
  else if (i .le. dim) then
    A(i,i) = (1.0 - 2.0*val)
  end if

  ! Fill top diagonal
  if ((i .gt. 1).and.(i .lt. (dim))) then
    A(i,i+1) = val
  end if
  ! Fill bottom diagonal
  if (i .lt. (dim-2)) then
    A(i+1,i) = val
  end if

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mv_mult(dim,A,V,O)
! Subroutine to multiply a square matrix and vector of corresponding size.

implicit none
  ! Dimension of matrix/vector
  integer,intent(in)   :: dim
  ! Define matrix and vector
  real,intent(in)      :: A(dim,dim), V(dim)
  ! Define output vector
  real,intent(out)     :: O(dim)

  ! Looping integers
  integer           :: i, j

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
