! Example code to demonstrate the finite differences schemes.
! forward, backwards, and central differences

! 9-3-19
! Jon Parsons

program finite_differences

implicit none
  integer   :: i ! looping integer
  real      :: h, a, b ! step, boundaries
  real      :: m, f_prime ! cells, derivative
  real      :: f, x ! Function variables, original
  real      :: p, y ! Function variables, derivative
  real      :: x1, x2, xc ! Temporary values

! Original
f(x) = cos(x)*log(x)
! Derivative (true)
p(y) = (cos(y)/y) - log(y)*sin(y)

! Initializations
m = 50.0
a = 1
b = 10

h = (b-a)/m

! true
open(unit=15,file="true.dat",status="replace",position="append")

tru_loop: do i = 1, int(m), 1
  x1 = h*float(i) + a

  write(15,*) x1, p(x1)

end do tru_loop

close(15)


! forward

open(unit=15,file="for.dat",status="replace",position="append")

for_loop: do i = 1, int(m), 1
  x1 = h*float(i) + a
  x2 = x1 + h

  f_prime = (f(x2) - f(x1))/h

  write(15,*) x1, f_prime

end do for_loop

close(15)

! backward

open(unit=15,file="back.dat",status="replace",position="append")

bac_loop: do i = 1, int(m), 1
  x1 = h*float(i) + a
  x2 = x1 + h

  f_prime = (f(x2) - f(x1))/h

  write(15,*) x2, f_prime

end do bac_loop

close(15)

! central

open(unit=15,file="cen.dat",status="replace",position="append")

cen_loop: do i = 2, int(m)-1, 1
  x1 = h*float(i-1) + a
  x2 = h*float(i+1) + a
  xc = h*float(i) + a

  f_prime = (f(x2) - f(x1))/(2*h)

  write(15,*) xc, f_prime

end do cen_loop

close(15)

end program
