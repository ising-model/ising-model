module montecarlo2d
    use hyperparameter
    use, intrinsic :: iso_fortran_env
    implicit none
    PUBLIC :: simulate2d
    PRIVATE

contains
    ! allocate memory to spin and initialize all spins to 1
    ! -----------------------------------------------------
    subroutine init_spin(spin)
        integer, dimension(:, :), allocatable :: spin
        allocate(spin(size, size))
        spin = 1
    end subroutine init_spin

    ! Metropolis-Hastings algorithm for Monte Carlo simulation
    ! --------------------------------------------------------
    subroutine metropolis2d(spin, beta)
        ! argument of the subroutine
        integer, dimension(:, :), allocatable :: spin
        real, intent(in)                      :: beta
        ! local variables of the subroutine
        integer i, x, y, s, R
        integer xpp, ypp, xnn, ynn ! neighbor spins
        real x0, y0, rand
        real dH

        ! Metropolis-Hastings algorithm main body
        do i = 1, size ** dim
            call random_number(x0)
            call random_number(y0)
            x = int(1 + (size - 1) * x0) ! x coordinate of selected spin
            y = int(1 + (size - 1) * y0) ! y coordinate of selected spin
            s = spin(x, y)               ! randomly selected spin

            ! obtain neighbor spins
            if (x == 1) then
                xpp = spin(2, y)
                xnn = spin(size, y)
            elseif (x == size) then
                xpp = spin(1, y)
                xnn = spin(size - 1, y)
            else
                xpp = spin(x + 1, y)
                xnn = spin(x - 1, y)
            endif

            if (y == 1) then
                ypp = spin(x, 2)
                ynn = spin(x, size)
            elseif (y == size) then
                ypp = spin(x, 1)
                ynn = spin(x, size - 1)
            else
                ypp = spin(x, y + 1)
                ynn = spin(x, y - 1)
            endif

            ! flip the spin and decide the next state
            R = xpp + ypp + xnn + ynn
            dH = 2 * s * R
            call random_number(rand)
            if (dH < 0.) then
                s = -s
            elseif (rand < exp(-beta * dH)) then
                s = -s
            endif
            spin(x, y) = s
        end do
    end subroutine metropolis2d

    ! calculation of energy of the spin using cshift
    ! ------------------------------
    real(real64) function calc_energy(spin)
        integer, dimension(:, :), allocatable, intent(in) :: spin
        integer, dimension(:, :), allocatable :: xps
        integer, dimension(:, :), allocatable :: yps
        integer, dimension(:, :), allocatable :: xns
        integer, dimension(:, :), allocatable :: yns

        ! shifted configuration of spins for easy computation
        xps = cshift(spin, shift=1, dim=1)
        yps = cshift(spin, shift=1, dim=2)
        xns = cshift(spin, shift=-1, dim=1)
        yns = cshift(spin, shift=-1, dim=2)

        calc_energy = sum(spin * -(xps + yps + xns + yns)) / (2 * dim)
    end function calc_energy

    ! calculation of magnetization of spin
    ! ------------------------------------
    real(real64) function calc_magnetization(spin)
        integer, dimension(:, :), allocatable, intent(in) :: spin
        calc_magnetization = sum(spin) / (2 * dim)
    end function calc_magnetization

    ! naive calculation of energy of the spin
    ! ---------------------------------
    real function calc_energy_naive(spin)
        integer, dimension(:, :), allocatable, intent(in) :: spin
        real(real64) En
        integer x, y, R
        integer xpp, ypp, xnn, ynn ! neighbor spins
        En = 0

        do x = 1, size
            do y = 1, size
                ! obtain neighbor spins
                if (x == 1) then
                    xpp = spin(2, y)
                    xnn = spin(size, y)
                elseif (x == size) then
                    xpp = spin(1, y)
                    xnn = spin(size - 1, y)
                else
                    xpp = spin(x + 1, y)
                    xnn = spin(x - 1, y)
                endif
    
                if (y == 1) then
                    ypp = spin(x, 2)
                    ynn = spin(x, size)
                elseif (y == size) then
                    ypp = spin(x, 1)
                    ynn = spin(x, size - 1)
                else
                    ypp = spin(x, y + 1)
                    ynn = spin(x, y - 1)
                endif

                R = xpp + ypp + xnn + ynn
                En = En - spin(x, y) * R
            end do
        end do
        calc_energy_naive = En / (2 * dim)
    end function calc_energy_naive

    subroutine simulate2d(beta, E, M, C, X, progress, step_progress)
        real, intent(in)    :: beta
        real, intent(inout) :: progress
        real, intent(in)    :: step_progress
        real, intent(inout) :: E, M, C, X
        real(real64)        :: tempE, tempM
        real(real64)        :: E1, M1, E2, M2
        real(real64)        :: norm, volume
        integer             :: step
        101 format(a, "Progress: ", f5.1, '% / ', f5.1, '%')

        ! initialize spin configuration
        integer, dimension(:, :), allocatable :: spin
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
            call metropolis2d(spin, beta)
            ! show progress bar
            !$OMP ATOMIC
            progress = progress + step_progress
            write (*, 101, advance='no') creturn, progress, 100.0
        end do
        ! Monte Carlo steps
        do step = 1, mcstep
            call metropolis2d(spin, beta)
            tempE = calc_energy(spin)
            tempM = calc_magnetization(spin)
            E1 = E1 + tempE / norm
            M1 = M1 + tempM / norm
            E2 = E2 + tempE ** 2 / norm
            M2 = M2 + tempM ** 2 / norm
            !$OMP ATOMIC
            progress = progress + step_progress
            write (*, 101, advance='no') creturn, progress, 100.0
        end do
        E = E1
        M = M1
        C = (E2 - volume * E1 ** 2) * beta ** 2
        X = (M2 - volume * M1 ** 2) * beta
    end subroutine simulate2d

end module