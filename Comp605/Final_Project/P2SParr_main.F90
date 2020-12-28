! Code written for the final project in COMP605.
! Adapted from the Phage-to-Shark project (Katira lab)

! Author: Jon Parsons
! 3-25-19

program P2SParr

! Declare modules
use mpi
use p2s_functions

implicit none

! model variables
	integer				:: gridsize, gridsize_local ! Original dimension; dimension of distributed arrays
	integer				:: tsteps, finished ! Number of iterations to run the simulation; flag for when
											! t = tsteps
	real,allocatable	:: coral(:,:), coral_temp(:,:) ! Original array and update array
	real,allocatable	:: coral_local(:,:), coral_temp_local(:,:) ! Arrays to distribute
	real				:: fish, shark, fish_cap ! fish population; predator population
												 ! fish carrying capacity
	real				:: corperc_ini, corperc  ! initial and current coral percentage of array

	integer				:: allocate_err ! Error variable for allocation
	integer				:: i, t ! Looping integers

! MPI variables
	integer				:: world_size, local_add ! number of processors; local address
	integer				:: mpi_err ! error variable

! Timing variables
	real(kind=8)		:: t1, t2 ! Holds time at start and time at end respectively

! Initialize MPI
call MPI_INIT (mpi_err)
call MPI_COMM_SIZE (MPI_COMM_WORLD,world_size, mpi_err)
call MPI_COMM_RANK (MPI_COMM_WORLD,local_add, mpi_err)

! Initialize flag
finished = 0

! Split tasks between processes
! Master process tasks
if (local_add .eq. 0) then

! Get inputs and initialize
	write(*,*) "Input number of timesteps: "
	read(*,*) tsteps
	write(*,*) "Input size of grid: "
	read(*,*) gridsize
	write(*,*) "Input initial coral percentage: "
	read(*,*) corperc_ini
	write(*,*) "Input initial fish pop: "
	read(*,*) fish
	write(*,*) "Input shark pop: "
	read(*,*) shark

	! Initialize local gridsize
	if (world_size .gt. 2) then
		gridsize_local = gridsize/world_size + 2
	else
		gridsize_local = gridsize
	end if

		! Initialize fish population
		fish_cap = fish

		! Initialize time
		t = 0

		! Allocation statements with initialization
		allocate(coral(gridsize,gridsize), stat=allocate_err)
			if (allocate_err .ne. 0) stop "Coral Allocation Failed"
		allocate(coral_temp(gridsize,gridsize), stat=allocate_err)
			if (allocate_err .ne. 0) stop "Coral_temp Allocation Failed"
		coral = 0.0
		coral_temp = 0.0

		! Populate initial coral layer
		write(*,*) "Populating coral layer"
		call initialize_coral(gridsize,coral,corperc_ini)
		corperc = percentcor(gridsize,coral)
		write(*,*) "Coral coverage at ", corperc, "%"

		! Distribute array sizes to daughter processes
		do i = 1, world_size-1, 1
			call MPI_SEND(gridsize,1,MPI_INT,i,0,MPI_COMM_WORLD,mpi_err)
			call MPI_SEND(gridsize_local,1,MPI_INT,i,0,MPI_COMM_WORLD,mpi_err)
		end do

		! Open file to hold timing data
		open(unit=17,file="timings.dat",status="unknown",position="append")
! Daughter proccess tasks
else

		! Recieve array dimensions
		call MPI_RECV(gridsize,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_err)
		call MPI_RECV(gridsize_local,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE, mpi_err)

		! Allocate arrays with initialization
		allocate(coral_local(gridsize,gridsize_local), stat=allocate_err)
			if (allocate_err .ne. 0) stop "Coral_local Allocation Failed"
		allocate(coral_temp_local(gridsize,gridsize_local), stat=allocate_err)
			if (allocate_err .ne. 0) stop "Coral_temp_local Allocation Failed"

		coral_local = 0.0
		coral_temp_local = 0.0


end if

! Get time of start
t1 = MPI_WTIME()


!!!!!!!!! Begin Iteration Loop !!!!!!!!!!!
! Ensure process synchronization
101	CALL MPI_Barrier(MPI_COMM_WORLD,mpi_err)

	! Distribute layer
	call coral_pass(gridsize,gridsize_local,coral,coral_local,local_add,world_size)

	! Update fish in master process
	if (local_add .eq. 0) then
		call pisc_update(fish,shark,fish_cap)
		write(*,*) "fish updated"
	end if

	! Distribute fish and carrying capacity variables
	call pass_real(fish,local_add,world_size)
	call pass_real(fish_cap,local_add,world_size)

	! Ensure synchronization S.T. each process has the requisite variable values
	CALL MPI_Barrier(MPI_COMM_WORLD,mpi_err)

	! Daughter processes update the coral layer
	if (local_add .ne. 0) then
		call coral_grow(gridsize,gridsize_local,coral_local)
		call coral_update(gridsize,gridsize_local,coral_local,coral_temp_local,fish,fish_cap)
	end if

	! Synchronize to keep processes aligned for sending local arrays back to master
	CALL MPI_Barrier(MPI_COMM_WORLD,mpi_err)

	! Return coral sub-arrays to master process for analysis
	call coral_recieve(gridsize,gridsize_local,coral_temp,coral_temp_local,local_add,world_size)

	! Master process updates time and displays current values
	if (local_add .eq. 0) then
		t = t + 1
		coral = coral_temp ! Set holding array to working array
		write(*,*) "Time", t
		write(*,*) "Coral Percentage: ", percentcor(gridsize,coral)
		write(*,*) "Fish pop: ", fish

		! Write current values to file
		call outputs(gridsize,coral,t,fish)

		! Check if the requested number of iterations has been met
		if (t .eq. tsteps) then
			finished = 1
			call MPI_BCAST(finished,1,MPI_INT,0,MPI_COMM_WORLD,mpi_err)
		end if
	end if

	! Check flag for value
	if (finished .eq. 1) then
		goto 103 ! Exit iteration loop
	else
		goto 101 ! Continue loop
	end if

!!!!!!!!!!!!!!!!!! End of Iteration Loop !!!!!!!!!!!!!!!!!!!
103 write(*,*) "Thread", local_add, "Has Finished instructions"

	! Get final time
    t2 = MPI_WTIME()

! Write time to completion to file
if (local_add .eq. 0) then
	write(17,*) t2-t1, world_size
end if


end program
