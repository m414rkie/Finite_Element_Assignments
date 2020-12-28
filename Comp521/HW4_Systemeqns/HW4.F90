program systemroot

implicit none
	real				:: temp(3), vars(3)
	real				:: jac(3,3), jacinv(3,3), root(3)
	real,allocatable	:: error(:,:)
	integer				:: maxiter, iter, i
	real				:: errtol


errtol = 0.0000001
maxiter = 1000
vars = 0.0
iter = 1
allocate(error(3,maxiter))

call funevals(vars,root)


open(unit=16,file="xyzerr.dat",status="replace",position="append")

do while ((abs(maxval(root)) .gt. errtol).and.(iter .lt. maxiter))

error(1,iter) = root(1)
error(2,iter) = root(2)
error(3,iter) = root(3)

call jacob(vars,jac)
call inverter(jac,jacinv,3)
call funevals(vars,root)


call multiplier(jacinv,root,temp)

vars = vars - temp
iter = iter + 1

end do

do i = 1, 3, 1
	error(i,:) = error(i,:) - root(i)
end do

do i = 1, iter, 1
	write(16,*) i, log(abs(error(1,i))), log(abs(error(2,i))), log(abs(error(3,i)))
end do

deallocate(error)
close(16)

end program


subroutine inverter(a,ainv,dime)
! Finds the inverse of a square matrix by solving a doubly-sized matrix
! as [A:I] => [I:A-1]

implicit none
	integer				:: dime
	real				:: a(dime,dime), ainv(dime,dime)
	integer				:: i,j, dub
	real,allocatable	:: work(:,:)

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

subroutine funevals(in,out)

implicit none
	real		:: in(3), out(3)
	real		::f1, f2, f3
	real		:: x, y, z

f1(x,y,z) = 3.0*x+5.0*y-z + 1.0/(x+2.0) - 1.0

f2(x,y,z) = 5.0*x+3.0*y-2.0*z + 1.0/(y+2.0) - 2.0

f3(x,y,z) = 3.0*x+5.0*y-3.0*z + 1.0/(z+2.0) - 3.0

out(1) = f1(in(1),in(2),in(3))
out(2) = f2(in(1),in(2),in(3))
out(3) = f3(in(1),in(2),in(3))

write(*,*) "X	Y	Z"
write(*,*) in(1), in(2), in(3)
write(*,*) "Function values"
write(*,*) out(1), out(2), out(3)

end subroutine

subroutine jacob(in,out)

implicit none
	real		:: in(3), out(3,3)
	real		:: f, x
	integer		:: i

f(x) = -1.0/((x+2.0)**2)

out(1,1) = 3.0
out(1,2) = 5.0
out(1,3) = 3.0
out(2,1) = 5.0
out(2,2) = 3.0
out(2,3) = 5.0
out(3,1) = -1.0
out(3,2) = -2.0
out(3,3) = -3.0

do i = 1, 3, 1
	out(i,i) = out(i,i) + f(in(i))
end do

end subroutine

subroutine multiplier(sqr,vec,out)

implicit none
	real	:: sqr(3,3)
	real	:: vec(3), out(3)
	integer :: i, j
out = 0.0
do i = 1, 3, 1

	do j = 1, 3, 1

		out(i) = out(i) + sqr(j,i)*vec(j)

	end do

end do

end subroutine
