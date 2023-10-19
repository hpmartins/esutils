subroutine readinput
use modinit
use modfunc
use modcore
implicit none

    integer iostat
    integer stmp(10)
    character(256) pname

    ! default values
    PMSTART=-1
    PMDELTA=-100.0
    PMDELTAC=0.0
    PMU=-100.0
    PMTS=-100.0
    PMTC=0.0
    PMDQ=-100.0
    PMJ=-100.0
    PMCOH=0
    PMEP=-5.0

    open(unit=80, file=PMINPFL, action='READ', iostat=iostat)
    if (iostat .ne. 0) then
        write(*,*)
        write(*,*) 'Unable to open input file: ', PMINPFL
        write(*,*)
        stop
    end if
10  continue
    read(80, *, end=30) pname
    if ((scan(trim(pname), '!') .eq. 1).or.(scan(trim(pname), '#') .eq. 1)) goto 10
    select case(trim(tolower(pname)))
    case('start')
        read(80, *, err=20) stmp(:)
        PMSTART = vec2bin(stmp)
    case('coherent')
        read(80, *, err=20) stmp(:)
        PMCOH = vec2bin(stmp)
    case('delta')
        read(80, *, err=20) PMDELTA
    case('u')
        read(80, *, err=20) PMU
    case('ts')
        read(80, *, err=20) PMTS
    case('dq')
        read(80, *, err=20) PMDQ
    case('j')
        read(80, *, err=20) PMJ
    case('ep')
        read(80, *, err=20) PMEP
    case('deltac')
        read(80, *, err=20) PMDELTAC
    case('tc')
        read(80, *, err=20) PMTC
    case('')
        goto 10
    end select
    goto 10
20  continue
    write(*,*)
    write(*,'("Problem occurred in ''",A,"'' pname")') trim(pname)
    write(*,*)
    stop
30  continue
    close(80)
    
    PMED = PMDELTA + PMEP - popcnt(PMSTART)*PMU
    PMEC = PMED - PMDELTAC + popcnt(PMSTART)*PMU
    PMTP = -0.5248*PMTS
    
    PMTLIST = [PMTS, PMTS, PMTP, PMTP, PMTP, PMTS, PMTS, PMTP, PMTP, PMTP]
end subroutine


subroutine basis
use modinit
use modfunc
use modcore
implicit none

    ! Ground state
    BGND = core_basis(PMSTART, PMCOH)
    gndsiz = size(BGND, 1)

    ! Addition states
    addlst = wbits(ieor(PMSTART, 1023))
    addsizup = 0
    addsizdn = 0
    do i = 1, size(addlst)
        j = addlst(i)
        call basissiz(siz, ibset(PMSTART, j), PMCOH)
        allocate(BADD(j) % v(siz, 5))
        BADD(j) % v = core_basis(ibset(PMSTART, j), PMCOH)
        
        if (j <= 4) then
            addsizdn = addsizdn + size(BADD(j)%v, 1)
        else
            addsizup = addsizup + size(BADD(j)%v, 1)
        end if
    end do

    ! Removal states
    remlst = wbits(PMSTART)
    remsizup = 0
    remsizdn = 0
    do i = 1, size(remlst)
        j = remlst(i)
        call basissiz(siz, ibclr(PMSTART, j), PMCOH)
        allocate(BREM(j) % v(siz, 5))
        BREM(j)%v = core_basis(ibclr(PMSTART, j), PMCOH)
        
        if (j <= 4) then
            remsizdn = remsizdn + size(BREM(j)%v, 1)
        else
            remsizup = remsizup + size(BREM(j)%v, 1)
        end if
    end do
end subroutine


subroutine hamiltonians
use modinit
use modfunc
use modcore
implicit none
    ! Ground state
    allocate(HGND(size(BGND, 1), size(BGND, 1))) 
    HGND = core_hamil(BGND)
    HCOR = HGND
    do i = 1, size(HCOR, 1)
        HCOR(i,i) = HCOR(i,i) - BGND(i,1)*PMU/0.83
    end do

    !Addition states
    do i = 1, size(addlst)
        j = addlst(i)
        allocate(HADD(j)%v(size(BADD(j)%v, 1), size(BADD(j)%v, 1)))
        HADD(j)%v = core_hamil(BADD(j)%v)
    end do

    !Removal states
    do i = 1, size(remlst)
        j = remlst(i)
        allocate(HREM(j)%v(size(BREM(j)%v, 1), size(BREM(j)%v, 1)))
        HREM(j)%v = core_hamil(BREM(j)%v)
    end do   
end subroutine


subroutine eigensolver
use modinit
use modfunc
use modcore
use types, only: dp
use linalg, only: eigh
implicit none
    ! Ground state
    real(dp), allocatable :: eva(:), evc(:,:)
    allocate(eva(size(HGND, 1)))
    allocate(evc(size(HGND, 1), size(HGND, 1)))
    call eigh(HGND, eva, evc)
    EVAGND = eva(1)
    EVCGND = evc(:, 1)
    
    ! Core state
    allocate(EVACOR(size(HCOR, 1)))
    allocate(EVCCOR(size(HCOR, 1), size(HCOR, 1)))
    call eigh(HCOR, EVACOR, EVCCOR)
    
    !Addition states
    do i = 1, size(addlst)
        j = addlst(i)
        allocate(EVAADD(j)%v(size(HADD(j)%v, 1)))
        allocate(EVCADD(j)%v(size(HADD(j)%v, 1), size(HADD(j)%v, 1)))
        call eigh(HADD(j)%v, EVAADD(j)%v, EVCADD(j)%v)
    end do
    
    !Removal states
    do i = 1, size(remlst)
        j = remlst(i)
        allocate(EVAREM(j)%v(size(HREM(j)%v, 1)))
        allocate(EVCREM(j)%v(size(HREM(j)%v, 1), size(HREM(j)%v, 1)))
        call eigh(HREM(j)%v, EVAREM(j)%v, EVCREM(j)%v)
    end do
end subroutine


subroutine operators
use modinit
use modfunc
use modcore
use types, only: dp
implicit none
    integer :: tvi(5), tvf(5), tvk(5)

    ! Addition states
    do i = 1, size(addlst)
        j = addlst(i)
        allocate(OPADDD(j)%v(size(BGND, 1), size(BADD(j)%v, 1)))
        allocate(OPADDP(j)%v(size(BGND, 1), size(BADD(j)%v, 1)))
        allocate(OPADDC(j)%v(size(BGND, 1), size(BADD(j)%v, 1)))
        OPADDD(j)%v(:,:) = 0.0_dp
        OPADDP(j)%v(:,:) = 0.0_dp
        OPADDC(j)%v(:,:) = 0.0_dp
        do b = 1, size(BADD(j)%v, 1)
            tvf = BADD(j)%v(b, :)
            do a = 1, size(BGND, 1)
                tvi = BGND(a, :)
                tvk = tvf - tvi
                if ((tvk(1) == 0) .and. (popcnt(tvk(3)) == 1) .and. &
                    (tvk(4) == 0) .and. (tvk(5) == 0)) then
                    OPADDP(j)%v(a,b) = 1.0_dp
                end if
                
                if ((tvk(1) == 1) .and. (tvk(3) == 0)) then
                    if ((tvk(4) == 0) .and. (tvk(5) == 0)) then
                        OPADDD(j)%v(a,b) = 1.0_dp
                    end if
                    if ((popcnt(tvk(4)) == 1) .and. (tvk(5) == 0)) then
                        OPADDC(j)%v(a,b) = 1.0_dp
                    end if
                    if ((tvk(4) == 0) .and. popcnt(tvk(5)) == 1) then
                        OPADDC(j)%v(a,b) = 1.0_dp
                    end if
                end if
            end do
        end do
    end do
    
    ! Removal states
    do i = 1, size(remlst)
        j = remlst(i)
        allocate(OPREMD(j)%v(size(BGND, 1), size(BREM(j)%v, 1)))
        allocate(OPREMP(j)%v(size(BGND, 1), size(BREM(j)%v, 1)))
        allocate(OPREMC(j)%v(size(BGND, 1), size(BREM(j)%v, 1)))
        OPREMD(j)%v(:,:) = 0.0_dp
        OPREMP(j)%v(:,:) = 0.0_dp
        OPREMC(j)%v(:,:) = 0.0_dp
        do b = 1, size(BREM(j)%v, 1)
            tvf = BREM(j)%v(b, :)
            do a = 1, size(BGND, 1)
                tvi = BGND(a, :)
                tvk = tvf - tvi
                if ((tvk(1) == 0) .and. (popcnt(tvk(3)) == 1) .and. &
                    (tvk(4) == 0) .and. (tvk(5) == 0)) then
                    OPREMP(j)%v(a,b) = 1.0_dp
                end if
                
                if ((tvk(1) == -1) .and. (tvk(3) == 0)) then
                    if ((tvk(4) == 0) .and. (tvk(5) == 0)) then
                        OPREMD(j)%v(a,b) = 1.0_dp
                    end if
                    if ((popcnt(tvk(4)) == 1) .and. (tvk(5) == 0)) then
                        OPREMC(j)%v(a,b) = 1.0_dp
                    end if
                    if ((tvk(4) == 0) .and. popcnt(tvk(5)) == 1) then
                        OPREMC(j)%v(a,b) = 1.0_dp
                    end if
                end if
            end do
        end do
    end do
end subroutine


subroutine transitions
use modinit
use modfunc
use modcore
use types, only: dp
implicit none

    integer               :: kup, kdn
    real(dp), allocatable :: tmrw(:,:)
    real(dp)              :: mxvl

    ! Core
    allocate(ACOR(gndsiz, 2))
    ACOR(:, 1) = EVACOR - minval(EVACOR)
    ACOR(:, 2) = (matmul(transpose(EVCCOR), EVCGND))**2
    
    ! Removal
    allocate(AREMUP(remsizup, 5))
    allocate(AREMDN(remsizdn, 5))
    AREMUP(:,:) = 0.0
    AREMDN(:,:) = 0.0
    kup = 1
    kdn = 1
    do i = 1, size(remlst)
        j = remlst(i)
        allocate(tmrw(size(BREM(j)%v,1), 5))
        tmrw(:, 1) = EVAGND - EVAREM(j)%v
        tmrw(:, 2) = (matmul(transpose(EVCREM(j)%v), matmul(transpose(OPREMD(j)%v), EVCGND)))**2
        tmrw(:, 3) = (matmul(transpose(EVCREM(j)%v), matmul(transpose(OPREMP(j)%v), EVCGND)))**2
        tmrw(:, 4) = (matmul(transpose(EVCREM(j)%v), matmul(transpose(OPREMC(j)%v), EVCGND)))**2
        tmrw(:, 5) = sum(tmrw(:, 2:4), dim=2)
        if (j >= 5) then
            AREMUP(kup:(kup+size(tmrw,1)-1), :) = tmrw
            kup = kup + size(tmrw,1)
        else
            AREMDN(kdn:(kdn+size(tmrw,1)-1), :) = tmrw
            kdn = kdn + size(tmrw,1)
        end if
        deallocate(tmrw)
    end do

    ! Addition
    allocate(AADDUP(addsizup, 5))
    allocate(AADDDN(addsizdn, 5))
    AADDUP(:,:) = 0.0
    AADDDN(:,:) = 0.0
    kup = 1
    kdn = 1
    do i = 1, size(addlst)
        j = addlst(i)
        allocate(tmrw(size(BADD(j)%v,1), 5))
        tmrw(:, 1) = EVAADD(j)%v - EVAGND
        tmrw(:, 2) = (matmul(transpose(EVCADD(j)%v), matmul(transpose(OPADDD(j)%v), EVCGND)))**2
        tmrw(:, 3) = (matmul(transpose(EVCADD(j)%v), matmul(transpose(OPADDP(j)%v), EVCGND)))**2
        tmrw(:, 4) = (matmul(transpose(EVCADD(j)%v), matmul(transpose(OPADDC(j)%v), EVCGND)))**2
        tmrw(:, 5) = sum(tmrw(:, 2:4), dim=2)
        if (j >= 5) then
            AADDUP(kup:(kup+size(tmrw,1)-1), :) = tmrw
            kup = kup + size(tmrw,1)
        else
            AADDDN(kdn:(kdn+size(tmrw,1)-1), :) = tmrw
            kdn = kdn + size(tmrw,1)
        end if
        deallocate(tmrw)
    end do
    
    
    
    mxvl = maxval([maxval(AREMUP(:,5)), maxval(AREMDN(:,5)), &
                   maxval(AADDUP(:,5)), maxval(AADDDN(:,5))])
    
    AREMUPT = filtertol(d_grpsum(AREMUP, 1), 5, 1e-3_dp*mxvl)
    AREMDNT = filtertol(d_grpsum(AREMDN, 1), 5, 1e-3_dp*mxvl)
    AADDUPT = filtertol(d_grpsum(AADDUP, 1), 5, 1e-3_dp*mxvl)
    AADDDNT = filtertol(d_grpsum(AADDDN, 1), 5, 1e-3_dp*mxvl)

end subroutine


subroutine export
use modinit
use modfunc
use modcore
use types, only: dp
implicit none

    open(unit=20, file=PMINPNM//'_ACOR.OUT', status='replace')
    do i = 1, size(ACOR, 1)
        write(20, *) ACOR(i, :)
    end do
    close(20)

    open(unit=20, file=PMINPNM//'_AREMUP.OUT', status='replace')
    do i = 1, size(AREMUPT, 1)
        write(20, *) AREMUPT(i, :)
    end do
    close(20)

    open(unit=20, file=PMINPNM//'_AREMDN.OUT', status='replace')
    do i = 1, size(AREMDNT, 1)
        write(20, *) AREMDNT(i, :)
    end do
    close(20)

    open(unit=20, file=PMINPNM//'_AADDUP.OUT', status='replace')
    do i = 1, size(AADDUPT, 1)
        write(20, *) AADDUPT(i, :)
    end do
    close(20)

    open(unit=20, file=PMINPNM//'_AADDDN.OUT', status='replace')
    do i = 1, size(AADDDNT, 1)
        write(20, *) AADDDNT(i, :)
    end do
    close(20)
end subroutine export