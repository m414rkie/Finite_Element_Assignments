! File containing subroutines for the parrallelized P2S project
! Contains subroutines on Lines;
! initialize_coral  : 13
! pass_real			: 37
! pass_int			: 68
! coral_pass		: 100
! coral_recieve		: 137
! pisc_update		: 171
! coral_grow		: 211
! coral_update		: 267
! outputs			: 381

subroutine initialize_coral(arrsize,cor_in,percentage)
! Subroutine to populate the initial coral layer

use p2s_functions

implicit none
	integer,intent(in)		:: arrsize ! Dimension of square input matrix
	real					:: cor_in(arrsize,arrsize) ! Input matrix
	real,intent(in)			:: percentage ! Desired percentage

! Array is filled with random numbers in range (0,1)
call random_number(cor_in)

! Set equal to zero all elements below the desired percentage range.
! This is not a perfect solution but works well with large arrays
where (cor_in .lt. percentage) cor_in = 0.0

! Scale remaining matrix to appropriate values
cor_in = cor_in*5.0

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pass_real(varpass,address,world)
! Subroutine to pass real variables from master process to
! daughter processes

! Module calls
use mpi

implicit none
	real	:: varpass ! Variable to pass
	integer :: address, world ! Local addresses; world size

	integer :: i, mpi_err ! Looping integer; error variable

if (address .eq. 0) then
! Master process distributes to all daughter processes

	do i = 1, world-1, 1
		call MPI_SEND(varpass,1,MPI_REAL,i,0,MPI_COMM_WORLD,mpi_err)
	end do

else
! Daughter process recieves sent variable

	call MPI_RECV(varpass,1,MPI_REAL,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_err)

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pass_int(varpass,address,world)
! Subroutine to pass integer variables from master process to
! daughter processes

! Module calls
use mpi

implicit none
	integer	:: varpass ! Variable to pass
	integer :: address, world ! Local addresses; world size

	integer :: i, mpi_err ! Looping integer; error variable

if (address .eq. 0) then
! Master process distributes to all daughter processes

	do i = 1, world-1, 1
		call MPI_SEND(varpass,1,MPI_INT,i,0,MPI_COMM_WORLD,mpi_err)
	end do

else
! Daughter processes recieve variable

	call MPI_RECV(varpass,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_err)

end if


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coral_pass(gridsize,gridsize_local,coral,coral_local,address,worldsize)
! Subroutine to break up and distribute the coral layer to daughter arrays

! Module call
use mpi

implicit none
	integer		:: gridsize, gridsize_local ! Initial dimension; local dimension
	real		:: coral(gridsize,gridsize), coral_local(gridsize,gridsize_local) ! Master grid; daughter grids

	integer		:: address, worldsize, mpi_err ! Local addresses; world size; error variable
	integer		:: i ! Looping integer

if (address .eq. 0) then
! Master process distributes portions of array to daughter processes
! Array is broken into column-delineated segments to take advantage of FORTRAN's column-major design

	do i = 1, worldsize-1, 1

		call MPI_SEND (coral(1:gridsize,(gridsize_local-1)*(i-1):(gridsize_local-1)*i),(gridsize*(gridsize_local)), &
						MPI_REAL,i,0,MPI_COMM_WORLD, mpi_err)
	end do

else
! Daughter processes recieve sent segment

	call MPI_RECV (coral_local,(gridsize*gridsize_local),MPI_REAL,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE, mpi_err)

	! Set lowest entries to zero. Assists with floating point error and removes excess coral
	where (coral_local .lt. 0.05) coral_local = 0.0

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coral_recieve(gridsize,gridsize_local,coral,coral_local,address,worldsize)
! Master process receives segments of main coral array from daughter processes

! Module calls
use mpi

implicit none
	integer		:: gridsize, gridsize_local ! Column and row dimension for local array
	real		:: coral(gridsize,gridsize), coral_local(gridsize,gridsize_local) ! working and holding arrays

	integer		:: address, worldsize, mpi_err ! Local address; world size; error variable
	integer		:: i ! Looping integer

! Daughter processes send local coral grid, tagged to ensure proper reassembly
do i = 1, worldsize-1, 1
	if (address .eq. i) then
			call MPI_SEND(coral_local(1:gridsize,2:(gridsize_local-1)),gridsize*(gridsize_local),MPI_REAL,0,i,MPI_COMM_WORLD,mpi_err)
	end if

end do

! Master process recieves sent data
if (address .eq. 0 ) then
	do i = 1, worldsize-1, 1
		call MPI_RECV(coral(1:gridsize,1+((i-1)*(gridsize_local-1)):(i*(gridsize_local-1))),gridsize*(gridsize_local), &
					MPI_REAL,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_err)
	end do

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pisc_update(prey,pred,carry_cap)
! Subroutine to update the fish.
! Update through simple logistic growth equation modified through predator interaction (stochastic)

implicit none
	real,intent(inout)		:: prey ! Input variable, fish
	real,intent(in)			:: carry_cap, pred ! Carrying capacity; predator population

	real					:: grow_rate, fish_delta ! Growth rate of prey; change in population
	real					:: pred_catch, day_avg, catch, eats ! Predator catch variable; average number of
																! days between sucessful predator hunts;
																! chance of success; amount removed from fish

! Initialize variables
grow_rate = 0.8
day_avg = 5.0

! Roll to determine successful hunt
call random_number(pred_catch)

! Set number to roll against
catch = 1.0/day_avg

! Hunt success logic
if (pred_catch .ge. catch) then
	eats = 4.0*pred/prey
else
	eats = 0.0
end if

! Determine amount of growth
fish_delta = grow_rate*(1.0 - (prey/carry_cap))*prey - pred*eats

! Update fish population
prey = prey + fish_delta

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coral_grow(gridsize,gridsize_local,coral_local)
! Subroutine relating to spawning or outward growth of coral

! Module call
use p2s_functions

implicit none
	integer,intent(in)		:: gridsize, gridsize_local ! column and row dimensions of array
	real,intent(inout)		:: coral_local(gridsize,gridsize_local) ! Local segment of coral

	integer					:: i, buds ! Looping integer; number of new corals
	real					:: coord(2) ! Array for holding coordinates of new coral
	integer					:: x, y ! Coordinates of new coral


! Determine number of new corals to place
buds = nint(sum(coral_local)/5.0)

! Coral can't grow that fast, limit growth
if (buds .gt. 30) then
	buds = 30
end if

! Working loop
do i = 1, buds, 1

! Find coordinates to place coral
102	call random_number(coord)

	! Scale to match grid
	x = floor(1.0 + coord(1)*gridsize)
	y = floor(1.0 + coord(2)*gridsize_local)

	! Ensure x, y are in bounds to avoid segfault
	if (x .gt. gridsize) then
		x = gridsize
	end if

	if (y .gt. gridsize_local) then
		y = gridsize_local
	end if

	! If no coral already present, place new coral
	if	(coral_local(x,y) .eq. 0.0) then
		coral_local(x,y) = 2.5
	else
		! If coral already present, return to coordinate call to try again
		goto 102
	end if

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine coral_update(gridsize,gridsize_local,coral_local,coral_temp_local,fish,fishcap)
! Subroutine to determine the growth or loss of the coral

implicit none
	integer,intent(in)		:: gridsize, gridsize_local ! Column and row dimensions of local array
	real,intent(in)			:: coral_local(gridsize,gridsize_local) ! Local coral array
	real,intent(in)			:: fish, fishcap ! Fish population; fish carrying capacity
	real					:: coral_temp_local(gridsize,gridsize_local) ! Holding array

	integer					:: i, j, widths, widthe ! Looping integers; endpoints of inner loop, handles low
												    ! processor counts

	real					:: growth_rate, alg_eats ! Natural coral growth rate; algae negative pressure on coral
	real					:: fish_impact ! Impact of fish on algae pressure
	integer					:: alg_count ! Holds count of algae around a particular point


! Logic to determine endpoints of row iterator
! At processor count = 2 the master array is passed in its entirety and padding is not needed
if (gridsize .eq. gridsize_local) then
	widths = 1
	widthe = gridsize
else
	widths = 2
	widthe = gridsize_local-1
end if

! Initialize variables
growth_rate = 1.08
alg_eats = 0.15/8.0
! Fish impact proportional to distance from carrying capacity
fish_impact = 0.1*(1.0 - (fish/fishcap))

! Ensure fish impact not negative and cannot inadvertantly assist algae
if (fish_impact .gt. 1.0) then
	fish_impact = 1.0
end if

if (fish_impact .lt. 0.0) then
	fish_impact = 0.0
end if

! Growth of coral
coral_temp_local = growth_rate*coral_local

! Ensure coral still within bounds of simulation ranges
where (coral_temp_local .gt. 5.0 ) coral_temp_local = 5.0

! Iterate through local coral for algae effect
do i = 1, gridsize, 1

	j_loop: do j = widths, widthe, 1

	! If not coral, skip this point
	if (coral_local(i,j) .eq. 0.0 ) then
		cycle j_loop
	end if

	! Initialize variable
	alg_count = 0

	! Logic statements ensure that the gridpoint in question has neighbors within bounds
	! and counts the algae around it. (counted algae is on holding array to ensure time consistency)
	if ((i .gt. 1) .and. (coral_local(i-1,j) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((i .lt. gridsize) .and. (coral_local(i+1,j) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((j .lt. gridsize_local) .and. (coral_local(i,j+1) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((j .gt. 1) .and. (coral_local(i,j-1) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((i .lt. gridsize) .and. (j .lt. gridsize_local) .and. (coral_local(i+1,j+1) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((i .gt. 1) .and. (j .gt. 1) .and. (coral_local(i-1,j-1) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((i .lt. gridsize) .and. (j .gt. 1) .and.(coral_local(i+1,j-1) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	if ((i .gt. 1) .and. (j .lt. gridsize_local) .and. (coral_local(i-1,j+1) .eq. 0.0)) then
		alg_count = alg_count + 1
	end if

	! If no algae around coral, negative effects are not felt, cycle loop
	if (alg_count .lt. 1) then
		cycle j_loop
	end if

	! Coral in question is lowered proportionally to the amount of algae and fish around it
	coral_temp_local(i,j) = coral_temp_local(i,j) - float(alg_count)*alg_eats*(1.0-fish_impact)

	end do j_loop

end do

! Remove coral below threshold
where (coral_temp_local .lt. 0.05) coral_temp_local = 0.0

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine outputs(gridsize,coral,time,fish)
! Subroutine to handle output at each timestep

! Module calls
use p2s_functions

implicit none
	integer,intent(in)	:: time, gridsize ! Current timestep; dimension of array
	real,intent(in)		:: coral(gridsize,gridsize), fish ! Input array; fish population

! Open files
open(unit=15,file="outs/coralperc.dat",status="unknown",position="append")
open(unit=16,file="outs/fishtime.dat",status="unknown",position="append")

! If timestep is first timestep, write headers to columns in files
if (time .eq. 1 ) then
	write(15,*) "Time percentage"
	write(16,*) "Time Fish"
end if

! Output to files
write(15,*) time, percentcor(gridsize,coral)
write(16,*) time, fish

end subroutine
