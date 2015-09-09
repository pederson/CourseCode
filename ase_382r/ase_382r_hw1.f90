program ase_382r_hw1

	implicit none

!	allocation part
	integer,parameter :: N = 1000000
	integer,parameter :: seed = 89758
	integer,parameter :: niters = 100
	integer,parameter :: Estart = 5;
	integer :: nparts = N/10;
	integer :: i, ncoll
	integer :: rval, p1, p2, hold
	integer :: coll_E_tot
	integer :: nt
	real :: S

	integer, dimension(:), allocatable :: coll_inds
	integer, dimension(:), allocatable :: En
	integer, dimension(:), allocatable :: f_E
	logical, dimension(:), allocatable :: will_collide

	allocate(coll_inds(nparts))
	allocate(En(N))
	allocate(f_E(N*Estart))
	allocate(will_collide(N))

	open(unit=5, file="entropy.txt")
	open(unit=1, file="dist_fn.txt")
	
	
	call srand(seed)

	! initialize the particle energies 
	forall (i=1:N) En(i) = Estart;

	nt = 1
	do while (nt<niters)

		! zero out the random subset of particles
		forall (i=1:N) will_collide(i) = .false.

		! select a random subset of particles that will collide (no repeats)
		ncoll = 0
		do while (ncoll < nparts)
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
		forall(i=1:N) f_E(i) = 0
		forall (i=1:N) f_E(En(i)) = f_E(En(i)) + 1

		!	print out the distribution function
		write (1, *) ' '
		write (1, "(A)", advance='no') 'TimeStep('
		write (1, "(I3)", advance='no') nt
		write (1, "(A)", advance='no') '): '
		i=1
		do while (i<50)
			write (1, "(I10)", advance='no') f_E(i)
			i = i+1
		enddo

		! compute the new entropy
		S = 0.0
		i=1
		do while (i .lt. N*Estart)
			if (f_E(i) .eq. 0) then
				i = i+1
				cycle
			endif
			if (i .gt. 500 .and. f_E(i) .eq. 0) exit
			S = S - real(f_E(i))/real(N) * log(real(f_E(i))/real(N))
			i=i+1
		enddo 

! 		print *, S

		! print out the entropy
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

end program ase_382r_hw1
