program ratiotest

implicit none
	integer :: i
	real	:: f, g, x, y
	real	:: ratiofit, ratiogit, stirling
	real	:: summf, summg, ratiof, ratiog
	
	
f(y) = 3.0**y

g(y) = 2.0*y/(2.0**y)

!ratiogit(y) = (y+1)/y

ratiofit = 3.0

i = 0

summf = f(float(i))
summg = g(float(i))
ratiof = f(float(i))
ratiog = g(float(i))

write(*,*) "Sum, Ratio, Sum, Ratio"
write(*,*) summf, ratiof, summg, ratiog

do i = 1, 15, 1

		summf = summf + f(float(i))
		summg = summg + g(float(i))
		ratiof = ratiof*3.0 + 1.0
		ratiog = ratiog + 1.0/(2.0**(i-1))
	
		write(*,*) summf, ratiof, summg, ratiog

end do

stirling = 50.0*log(50.0) - 50.0 + log(sqrt(2*3.1459*50.0))
write(*,*) stirling

end program