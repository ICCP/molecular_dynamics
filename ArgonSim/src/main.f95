program ArgonSim

    use initial_conditions
    use verlet
    !use sub_regions
    implicit none

    ! System parameters
    integer(16), parameter :: N = 500! 4, 32, 108, 256, 500, 864, 1372, 2048, 2916, 4000
    real(16), parameter :: SIGMA = 1
    real(16), parameter :: EPS = 1
    real(16), parameter :: MASS = 1
    real(16), parameter :: Kb = 1
    real(16), parameter :: DT = .004
    integer(8), parameter :: MAXSTEPS = 2000
    integer(8), parameter :: WRITEFREQ = 4
    integer(8), parameter :: TEMPFREQ = 20

    ! System array structure
    real(16) :: T = 0
    real(16), dimension(3, N) :: r, v, a, olda
    real(16), dimension(MAXSTEPS) :: temps
    real(16), dimension(4) :: temps_array
    real(16), dimension(200) :: r_histo
    real(16) :: one_cube_length, l, Ek, lambda, currTemp, vsquared
    real(16) :: Ep, mom(3), pressure, temp_a(3), dr
    integer(8) :: i, frame, k, cube_side, u
    integer(8) :: thermostat_bool = 0
    character(5) :: strframe

    temps_array(1) = 0
    temps_array(2) = 1
    temps_array(3) = 10

    open(unit = 1001, file = "data/energy.dat")
    open(unit = 1002, file = "data/temp.dat")
    open(unit = 1003, file = "data/pressure.dat")
    open(unit = 1004, file = "data/pair_correlation.dat")

do u = 1, 3
    T = temps_array(u)
    currTemp = T
    frame = 0 ! Counter frame used to generate positions data
    one_cube_length = 2d0**(2d0/3) * SIGMA ! Determine the volume of the system and its size
    cube_side = nint(((N*2d0)**(1d0/3d0))/2d0) ! Number of cubes on one side
    l = 1.0d0*one_cube_length*cube_side ! Length of side of available volume
    dr = SIGMA/10d0 ! Size of delta r used for pair correlation calculations
    !num_bins = nint(sqrt(3*l**2d0)/dr) + 1
    if (N /= ((2*cube_side)**(3))/2) then
        print *, "*****************"
        print *, "WARNING: Invalid number of particles."
        print *, "*****************"
    end if

    !allocate(r_histo(num_bins))
    r_histo = 0
    ! Initialize random distributions & lattices
    call max_boltz(MASS, T, v)
    call fcc_lattice(SIGMA, r, cube_side)

    ! Calculate initial energies and initialize accelerations
    call calc_energy(v, MASS, Ek, Ep, mom, SIGMA, EPS, r, l) ! From verlet_algorithm.f95
    call total_force(EPS, SIGMA, MASS, l, r, a, Ep, kB, currTemp, pressure, temp_a) ! From verlet_algorithm.f95

    do i = 1, MAXSTEPS ! Begin main loop
        !f( T == 5) then
            ! Monitor the results printing positions and energies to files
            if (mod(i,WRITEFREQ) == 1) then
                print *, "Timestep ", i, " of ",MAXSTEPS
                frame = frame + 1
                !print *, frame, "Ek: ", Ek, "Ep: ", Ep, "E: ", Ek+Ep
                write(strframe,'(I5)') frame
                open(unit = 1000, file = "data/position"//strframe//".dat")
                do k = 1, N
                  write(1000, *) r(:,k)
                end do
                close(1000)
            end if
        !nd if

        ! Runs the verlet time evolution for each time step
        call verlet_algorithm(EPS, SIGMA, MASS, l, dt, r, v, a, Ep, vsquared, olda, temp_a, pressure, currTemp, kB) ! From verlet_algorithm.f95
        call calc_energy(v, MASS, Ek, Ep, mom, SIGMA, EPS, r, l) ! From verlet_algorithm.f95
        

        ! THERMOSTAT
        currTemp = 2d0*Ek/(3d0*Kb)
        if (Ek < -10*Ep) then
            write(1001, *) Ek
            write(1001, *) Ep
            write(1001, *) Ek+Ep
        end if
        
        write(1003, *) pressure

        if (mod(i,TEMPFREQ) == 1) then
            !print *, "Desired temperature: ", T
            !print *, "Current temperature: ", currTemp
            lambda = (3d0*Kb*T/(MASS*vsquared))**0.5
            v=v*lambda/(N**.5)
            call calc_energy(v, MASS, Ek, Ep, mom, SIGMA, EPS, r, l)
            currTemp = 2d0*Ek/(3d0*Kb)
            !print *, "New temperature: ", currTemp
        end if

        if (currTemp < 10*T) then
            write(1002, *) currTemp
        end if

    end do

    call pair_correlation(dr, r, r_histo, l)

    do i = 1, size(r_histo)
        write(1004, *) r_histo(i)/((N-1)*N/2)
    end do

    write(1000, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    write(1001, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    write(1002, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    write(1003, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    write(1004, *) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
end do
    close(1001)
    close(1002)
    close(1003)
    close(1004)
end program ArgonSim

