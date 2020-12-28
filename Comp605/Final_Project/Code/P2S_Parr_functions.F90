module p2s_functions

contains

real function percentcor(size,arrin)
! Function takes square array as input and 
! returns the percentage of that array which
! is not zero

implicit none
	integer		:: size ! Dimension of input matrix
	real		:: arrin(size,size) ! Input matrix
	integer		:: corcount ! Holds a count of non-zero elements
	logical		:: wherecor(size,size) ! Logical array for use in finding zeroes
	
! Assign .TRUE. to elements corresponding to non-zero elements
wherecor = (arrin .gt. 0.0)

! Count .TRUE. elements
corcount = count(wherecor)

! Return percentage
percentcor = float(corcount)/(float(size)**2)

end function

end module