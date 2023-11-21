program main
    use OMP_LIB
    use f90getopt
    use hyperparameter
    use montecarlo3d
    implicit none

    ! hyperparameters entered as arguments
    ! ------------------------------------
    type(option_s) :: opts(8)

    ! statistical values calculated through Monte Carlo
    !--------------------------------------------------
    real, dimension(:), allocatable :: Ts, Es, Ms, Cs, Xs
    integer                         :: num_steps, index
    real                            :: T, E, M, C, X

    opts(1) = option_s('size', .true.,  's')
    opts(2) = option_s('dim',  .true., 'd')
    opts(3) = option_s('init_temp', .true., 'i')
    opts(4) = option_s('final_temp', .true., 'f')
    opts(5) = option_s('temp_step', .true., 't')
    opts(6) = option_s('mcstep', .true., 'm')
    opts(7) = option_s('eqstep', .true., 'e')
    opts(8) = option_s('help', .false., 'h')

    do 
        select case (getopt("s:d:i:f:t:m:e:h", opts))
            case ('s')
                read (optarg, '(i4)') size
            case ('d')
                read (optarg, '(i1)') dim
            case ('i')
                read (optarg, '(f7.5)') init_temp
            case ('f')
                read (optarg, '(f7.5)') final_temp
            case ('t')
                read (optarg, '(f7.6)') temp_step
            case ('m')
                read (optarg, '(i7)') mcstep
            case ('e')
                read (optarg, '(i7)') eqstep
            case ('h')
                call print_help()
                stop
            case default
                exit
        end select
    end do

    ! allocate array to store simulation results
    ! ------------------------------------------
    num_steps = int((final_temp - init_temp) / temp_step) + 1
    allocate(Ts(num_steps))
    allocate(Es(num_steps))
    allocate(Ms(num_steps))
    allocate(Cs(num_steps))
    allocate(Xs(num_steps))
    Ts = 0
    Es = 0
    Ms = 0
    Cs = 0
    Xs = 0

    ! parallel Monte Carlo simulation
    ! -------------------------------
    !$OMP PARALLEL DO PRIVATE(index, T, E, M, C, X)
    do index = 1, num_steps
        T = init_temp + temp_step * (index - 1)
        E = index
        M = 0
        C = 0
        X = 0

        call simulate3d(1.0 / T, E, M, C, X)
        Ts(index) = T
        Es(index) = E
        Ms(index) = M
        Cs(index) = C
        Xs(index) = X
    end do
    !$OMP END PARALLEL DO

    deallocate(Ts)
    deallocate(Es)
    deallocate(Ms)
    deallocate(Cs)
    deallocate(Xs)

contains
    subroutine print_help()
        print '(a, /)', 'command-line options:'
        print '(a)', '  -s, --size        size of the lattice               (default: 30)'
        print '(a)', '  -d, --dim         dimension of the lattice          (default: 3)'
        print '(a)', '  -i, --init_temp   initial temperature of the output (default: 1.5)'
        print '(a)', '  -f, --final_temp  final temperature of the output   (default: 6.5)'
        print '(a)', '  -t, --temp_step   step size of the temperature      (default: 0.04)'
        print '(a)', '  -m, --mcstep      number of Monte Carlo steps       (default: 1000)'
        print '(a)', '  -e, --eqstep      number of steps for equilibration (default: 1000)'
        print '(a, /)', '  -h, --help        print usage information and exit'
    end subroutine print_help
end program main