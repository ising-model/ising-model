module montecarlo3d
    use hyperparameter
    use, intrinsic :: iso_fortran_env
    implicit none
    PUBLIC :: simulate3d
    PRIVATE

contains
    ! allocate memory to spin and initialize all spins to 1
    ! -----------------------------------------------------
    subroutine init_spin(spin)
        integer, dimension(:, :, :), allocatable :: spin
        real, dimension(:, :, :), allocatable :: rand
        allocate(spin(size, size, size))
        allocate(rand(size, size, size))
        call random_number(rand)
        spin = 2 * int(2 * rand) - 1
    end subroutine init_spin

    ! Metropolis-Hastings algorithm for Monte Carlo simulation
    ! --------------------------------------------------------
    subroutine metropolis3d(spin, beta)
        ! argument of the subroutine
        integer, dimension(:, :, :), allocatable :: spin
        real, intent(in)                         :: beta
        ! local variables of the subroutine
        integer i, x, y, z, s, R
        integer xpp, ypp, zpp, xnn, ynn, znn ! neighbor spins
        real x0, y0, z0, rand
        real dH

        ! Metropolis-Hastings algorithm main body
        do i = 1, size ** dim
            call random_number(x0)
            call random_number(y0)
            call random_number(z0)
            x = int(1 + (size - 1) * x0) ! x coordinate of selected spin
            y = int(1 + (size - 1) * y0) ! y coordinate of selected spin
            z = int(1 + (size - 1) * z0) ! z coordinate of selected spin
            s = spin(x, y, z)            ! randomly selected spin

            ! obtain neighbor spins
            if (x == 1) then
                xpp = spin(2, y, z)
                xnn = spin(size, y, z)
            elseif (x == size) then
                xpp = spin(1, y, z)
                xnn = spin(size - 1, y, z)
            else
                xpp = spin(x + 1, y, z)
                xnn = spin(x - 1, y, z)
            endif

            if (y == 1) then
                ypp = spin(x, 2, z)
                ynn = spin(x, size, z)
            elseif (y == size) then
                ypp = spin(x, 1, z)
                ynn = spin(x, size - 1, z)
            else
                ypp = spin(x, y + 1, z)
                ynn = spin(x, y - 1, z)
            endif

            if (z == 1) then
                zpp = spin(x, y, 2)
                znn = spin(x, y, size)
            elseif (z == size) then
                zpp = spin(x, y, 1)
                znn = spin(x, y, size - 1)
            else
                zpp = spin(x, y, z + 1)
                znn = spin(x, y, z - 1)
            endif

            ! flip the spin and decide the next state
            R = xpp + ypp + zpp + xnn + ynn + znn
            dH = 2 * s * R
            call random_number(rand)
            if (dH < 0.) then
                s = -s
            elseif (rand < exp(-beta * dH)) then
                s = -s
            endif
            spin(x, y, z) = s
        end do
    end subroutine metropolis3d

    ! calculation of energy of the spin using cshift
    ! ------------------------------
    real function calc_energy(spin)
        integer, dimension(:, :, :), allocatable, intent(in) :: spin
        integer, dimension(:, :, :), allocatable :: xps
        integer, dimension(:, :, :), allocatable :: yps
        integer, dimension(:, :, :), allocatable :: zps
        integer, dimension(:, :, :), allocatable :: xns
        integer, dimension(:, :, :), allocatable :: yns
        integer, dimension(:, :, :), allocatable :: zns

        ! shifted configuration of spins for easy computation
        xps = cshift(spin, shift=1, dim=1)
        yps = cshift(spin, shift=1, dim=2)
        zps = cshift(spin, shift=1, dim=3)
        xns = cshift(spin, shift=-1, dim=1)
        yns = cshift(spin, shift=-1, dim=2)
        zns = cshift(spin, shift=-1, dim=3)

        calc_energy = sum(spin * -(xps + yps + zps + xns + yns + zns)) / (2 * dim)
    end function calc_energy

    ! calculation of magnetization of spin
    ! ------------------------------------
    real function calc_magnetization(spin)
        integer, dimension(:, :, :), allocatable, intent(in) :: spin
        calc_magnetization = sum(spin)
    end function calc_magnetization

    ! naive calculation of energy of the spin
    ! ---------------------------------
    real function calc_energy_naive(spin)
        integer, dimension(:, :, :), allocatable, intent(in) :: spin
        real En
        integer x, y, z, R
        integer xpp, ypp, zpp, xnn, ynn, znn ! neighbor spins
        En = 0

        do x = 1, size
            do y = 1, size
                do z = 1, size
                    ! obtain neighbor spins
                    if (x == 1) then
                        xpp = spin(2, y, z)
                        xnn = spin(size, y, z)
                    elseif (x == size) then
                        xpp = spin(1, y, z)
                        xnn = spin(size - 1, y, z)
                    else
                        xpp = spin(x + 1, y, z)
                        xnn = spin(x - 1, y, z)
                    endif

                    if (y == 1) then
                        ypp = spin(x, 2, z)
                        ynn = spin(x, size, z)
                    elseif (y == size) then
                        ypp = spin(x, 1, z)
                        ynn = spin(x, size - 1, z)
                    else
                        ypp = spin(x, y + 1, z)
                        ynn = spin(x, y - 1, z)
                    endif

                    if (z == 1) then
                        zpp = spin(x, y, 2)
                        znn = spin(x, y, size)
                    elseif (z == size) then
                        zpp = spin(x, y, 1)
                        znn = spin(x, y, size - 1)
                    else
                        zpp = spin(x, y, z + 1)
                        znn = spin(x, y, z - 1)
                    endif

                    R = xpp + ypp + zpp + xnn + ynn + znn
                    En = En - spin(x, y, z) * R
                end do
            end do
        end do
        calc_energy_naive = En / (2 * dim)
    end function calc_energy_naive

    subroutine simulate3d(beta, E, M, C, X)
        real, intent(in)    :: beta
        real, intent(inout) :: E, M, C, X
        real                :: tempE, tempM
        real(real64)        :: E1, M1, E2, M2
        integer             :: norm, volume, step

        ! initialize spin configuration
        integer, dimension(:, :, :), allocatable :: spin
        call init_spin(spin)
        E1 = 0
        M1 = 0
        E2 = 0
        M2 = 0

        ! Monte Carlo simulation
        volume = size ** dim
        norm = mcstep * volume
        ! equilibration steps
        do step = 1, eqstep
            call metropolis3d(spin, beta)
        end do
        ! Monte Carlo steps
        do step = 1, mcstep
            call metropolis3d(spin, beta)
            tempE = calc_energy(spin)
            tempM = calc_magnetization(spin)
            E1 = E1 + tempE
            M1 = M1 + tempM
            E2 = E2 + tempE ** 2 / norm
            M2 = M2 + tempM ** 2 / norm
        end do
        E = E1 / norm
        M = M1 / norm
        C = (E2 - volume * (E1 / norm) ** 2) * beta ** 2
        X = (M2 - volume * (M1 / norm) ** 2) * beta
    end subroutine simulate3d

end module