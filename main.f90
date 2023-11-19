program main
    implicit none
    character(len=128) :: arg
    integer            :: i

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)

        select case (arg)
            case ('-i', '--init_temp')
                print '(2a)', '1.0'
        end select
    end do
    
end program main