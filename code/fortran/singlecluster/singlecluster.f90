program singlecluster
use modinit
use modfunc
use modcore
implicit none
    narg = command_argument_count()
    call get_command_argument(0, EXCNM)
    if (narg .eq. 0) then
        inquire(file='singlecluster.in', exist=bool)
        if (bool) then
            PMINPFL = 'singlecluster.in'
            PMINPNM = 'singlecluster'
        else
            write(*,*)
            write(*,*) 'No input file given. Usage: ', trim(EXCNM), ' <inputfile>'
            write(*,*)
            stop
        end if
    else
        call get_command_argument(1, tmpstring)
        PMINPFL = trim2(tmpstring)
        PMINPNM = PMINPFL(1:(scan(PMINPFL, '.', BACK=.true.)-1))
    end if
    
    write(*,*) 'Input'
    call readinput
    write(*,*) 'Basis'
    call basis
    write(*,*) 'Hamiltonians'
    call hamiltonians
    write(*,*) 'Transition operators'
    call operators
    write(*,*) 'Eigenvalues and eigenvectors'
    call eigensolver
    write(*,*) 'Transitions'
    call transitions
    write(*,*) 'Exporting'
    call export
    write(*,*) 'Ok!'
    read *
end program singlecluster
