program main
    use OMP_LIB
    use f90getopt
    use hyperparameter
    use montecarlo2d
    use montecarlo3d
    implicit none

    ! hyperparameters entered as arguments
    ! ------------------------------------
    type(option_s) :: opts(9)

    ! statistical values calculated through Monte Carlo
    !--------------------------------------------------
    real, dimension(:), allocatable :: Ts, Es, Ms, Cs, Xs
    integer                         :: num_steps, index
    real                            :: progress
    real                            :: T, E, M, C, X
    character(len=7)                :: size_char, dim_char, eqstep_char, mcstep_char
    character(len=100)              :: csv_file

    ! receive argument from the command line
    ! --------------------------------------
    opts(1) = option_s('size', .true.,  's')
    opts(2) = option_s('dim',  .true., 'd')
    opts(3) = option_s('init_temp', .true., 'i')
    opts(4) = option_s('final_temp', .true., 'f')
    opts(5) = option_s('temp_step', .true., 't')
    opts(6) = option_s('mcstep', .true., 'm')
    opts(7) = option_s('eqstep', .true., 'e')
    opts(8) = option_s('dir', .true., 'r')
    opts(9) = option_s('help', .false., 'h')

    do 
        select case (getopt("s:d:i:f:t:m:e:r:h", opts))
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
            case ('r')
                read (optarg, '(a)') directory
            case ('h')
                call print_help()
                stop
            case default
                exit
        end select
    end do

    ! print number of threads for simulation
    ! --------------------------------------
    call print_num_threads()

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
    progress = 0.
    write (*, 101, advance='no') creturn, progress, 100.0
    101 format(a, "Progress: ", f5.1, '% / ', f5.1, '%')
    !$OMP PARALLEL DO PRIVATE(index, T, E, M, C, X)
    do index = 1, num_steps
        T = init_temp + temp_step * (index - 1)

        if (dim == 2) then
            call simulate2d(1.0 / T, E, M, C, X)
        else
            call simulate3d(1.0 / T, E, M, C, X)
        endif

        Ts(index) = T
        Es(index) = E
        Ms(index) = abs(M)
        Cs(index) = C
        Xs(index) = X

        ! show progress bar
        !$OMP ATOMIC
        progress = progress + 100. / num_steps
        write (*, 101, advance='no') creturn, progress, 100.0
    end do
    !$OMP END PARALLEL DO
    print *, ""

    ! save results to the directory
    ! -----------------------------
    write (size_char, '(i4.4)') size
    write (dim_char, '(i1)') dim
    write (eqstep_char, '(i7.7)') eqstep
    write (mcstep_char, '(i7.7)') mcstep
    csv_file = trim(directory) // trim("result_L") // trim(size_char) &
                               // trim("_D") // trim(dim_char) &
                               // trim("_EQ") // trim(eqstep_char) &
                               // trim("_MC") // trim(mcstep_char) // trim(".csv")
    open(unit=1, file=csv_file, status='unknown')
    write (1, '(5(a, x))') 'T', 'E', 'M', 'C', 'X'
    do index = 1, num_steps
        write (1, '(5(g20.10, x))') Ts(index), Es(index), Ms(index), Cs(index), Xs(index)
    end do
    close(1)
    print '(2a, /)', "Results saved to ", csv_file

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
        print '(a)', '  -r, --dir         directory to save the results     (default: ./results/)'
        print '(a, /)', '  -h, --help        print usage information and exit'
    end subroutine print_help

    subroutine print_num_threads()
        integer num_threads
        !$OMP PARALLEL SHARED(num_threads)
        num_threads = omp_get_num_threads()
        !$OMP END PARALLEL
        print '(a, i3)', "Number of threads: ", num_threads
    end subroutine print_num_threads
end program main