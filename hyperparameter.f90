module hyperparameter
    implicit none
    save
    integer           :: size = 30
    integer           :: dim = 3
    real              :: init_temp = 1.5
    real              :: final_temp = 6.5
    real              :: temp_step = 0.04
    integer           :: mcstep = 1000
    integer           :: eqstep = 1000
    character(len=20) :: directory = "./results/"
    character*1       :: creturn = achar(13)
end module