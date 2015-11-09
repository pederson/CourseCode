program ase_382r_hw1

	implicit none

!	allocation part
	integer,parameter :: N = 1000000
	integer,parameter :: seed = 89758
	real,parameter :: S_limit = 0.0001
	integer,parameter :: Estart = 5; ! 20
	integer :: nparts = N/10;
	integer :: i, ncoll
	integer :: rval, p1, p2, hold
	integer :: coll_E_tot
	integer :: nt
	real :: S, delta_S, S_old

	integer, dimension(:), allocatable :: coll_inds
	integer, dimension(:), allocatable :: En
	real, dimension(:), allocatable :: f_E
	logical, dimension(:), allocatable :: will_collide

	allocate(coll_inds(nparts))
	allocate(En(N))
	allocate(f_E(N*Estart))
	allocate(will_collide(N))

	open(unit=5, file="entropy.txt")
	open(unit=1, file="dist_fn.txt")
	
	
	call srand(seed)

	! initialize the particle energies 
 	forall (i=1:N) En(i) = Estart; ! single initial energy
	
 	! spread of initial energy
! 	forall(i=1:N/5) En(i) = Estart-2;
! 	forall(i=N/5:2*N/5) En(i) = Estart-1;
! 	forall(i=2*N/5:3*N/5) En(i) = Estart;
! 	forall(i=3*N/5:4*N/5) En(i) = Estart+1;
! 	forall(i=4*N/5:N) En(i) = Estart+2;

	! initialize the entropy and delta
	S = 0.1
	delta_S = 1.0
	S_old = S

	nt = 1
	do while (delta_S/S .gt. S_limit)

		! zero out the random subset of particles
		forall (i=1:N) will_collide(i) = .false.

		! select a random subset of particles that will collide (no repeats)
		ncoll = 1
		do while (ncoll .le. nparts)
			rval = rand()*(N-1) + 1;
			if (.not. will_collide(rval)) then
				will_collide(rval) = .true.
				coll_inds(ncoll) = rval
				ncoll = ncoll + 1;
			endif
		enddo


		ncoll = 0;
		do while (ncoll < nparts)
			! select random pairs from the subset
			if (nparts-ncoll<3) exit
			! inds between 1 and nparts-ncoll
			p1 = rand()*(nparts-ncoll-1) + 1
			p2 = rand()*(nparts-ncoll-1) + 1
			if (p1 .eq. p2) cycle
			if (p1 .eq. 0 .or. p2 .eq. 0) print *, 'zero index!'
			if (p1 .gt. N .or. p2 .gt. N) print *, 'seg fault'

			! perform the collision
			! compute total energy of the pairs
			coll_E_tot = En(coll_inds(p1)) + En(coll_inds(p2))
			! first particle leaves with a random fraction of the total energy
			En(coll_inds(p1)) = rand()*coll_E_tot
			En(coll_inds(p2)) = coll_E_tot - En(coll_inds(p1))

			! swap the two indices to the back of the index array
			hold = coll_inds(nparts-ncoll)
			coll_inds(nparts-ncoll) = coll_inds(p1)
			coll_inds(p1) = hold

			hold = coll_inds(nparts-ncoll-1)
			coll_inds(nparts-ncoll-1) = coll_inds(p2)
			coll_inds(p2) = hold

			ncoll = ncoll + 2
		enddo

		! compute the new energy distribution
		forall(i=1:N*Estart) f_E(i) = 0
		forall (i=1:N) f_E(En(i)+1) = f_E(En(i)+1) + 1
		forall(i=1:N*Estart) f_E(i) = f_E(i)/real(N)

		!	print out the distribution function to file
		write (1, *) ' '
		write (1, "(A)", advance='no') 'TimeStep('
		write (1, "(I3)", advance='no') nt
		write (1, "(A)", advance='no') '): '
		i=1
		do while (i<50)
			write (1, "(F11.6)", advance='no') f_E(i)
			i = i+1
		enddo

		! compute the new entropy
		S_old = S
		S = 0.0
		i=1
		do while (i .lt. N*Estart)
			if (f_E(i) .eq. 0.0) then ! dangerous comparison to float
				i = i+1
				cycle
			endif
			!if (i .gt. 5000 .and. f_E(i) .eq. 0) exit
			S = S - f_E(i) * log(f_E(i))
			i=i+1
		enddo 

		delta_S = S - S_old

! 		print *, S

		! print out the entropy to file
		write (5, "(A)", advance='no') 'TimeStep('
		write (5, "(I3)", advance='no') nt
		write (5, "(A)", advance='no') '): '
		write (5, *) S

		nt = nt + 1
	enddo



	! close files
	close(1, status='keep')
	close(5, status='keep')

	print *, 'THE END'

	deallocate(coll_inds)
 	deallocate(En)
	deallocate(f_E)
	deallocate(will_collide)

end program ase_382r_hw1
