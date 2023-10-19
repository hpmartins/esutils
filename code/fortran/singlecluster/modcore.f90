module modcore
use modinit
use modfunc
implicit none

contains

    function core_basis(st, coh) result(fout)
        integer, intent(in) :: st, coh
        integer, allocatable :: fout(:,:)
        integer i, j, k, m, nst
        integer coh1cmp, coh1cnt, defcnt, coh2cmp, coh2cnt, totcnt
        integer, allocatable :: coh1pos(:), coh2pos(:), buf(:)
        
        nst = popcnt(st)
        
        ! Total of configurations
        defcnt  = 2**(10 - nst)
        coh1cmp = iand(st, coh)
        coh1cnt = popcnt(coh1cmp)
        coh2cmp = iand(not(st), coh)
        coh2cnt = popcnt(coh2cmp)*(2**(10 - nst - 1))
        coh2pos = wbits(coh2cmp)
        totcnt  = defcnt + coh1cnt + coh2cnt
        allocate(fout(totcnt, 5))
        
        j = 1
        do i = 1, 1024
            if ((iand(i, st) == st)) then
                fout(j,1) = popcnt(i)
                fout(j,2) = i
                fout(j,3) = iand(i, not(st))
                fout(j,4) = 0
                fout(j,5) = 0
                j = j + 1
                
                if (coh2cnt .gt. 0) then
                    do m = 1, size(coh2pos)
                        k = coh2pos(m)
                        if (btest(i, k) .and. (.not. btest(st, k))) then
                            fout(j,1) = popcnt(i)
                            fout(j,2) = i
                            fout(j,3) = iand(i, not(ibset(st, k)))
                            fout(j,4) = 2**k
                            fout(j,5) = 0
                            j = j + 1
                        end if
                    end do
                end if
             end if
        end do

        if (coh1cnt .gt. 0) then
            coh1pos = wbits(coh1cmp)
            do i = 1, size(coh1pos)
                k = coh1pos(i)
                fout(j,1) = popcnt(ibclr(st, k))
                fout(j,2) = ibclr(st, k)
                fout(j,3) = 0
                fout(j,4) = 0
                fout(j,5) = 2**k
                j = j + 1
            end do
        end if

        allocate(buf(size(fout, 2)))
        do i = 1, size(fout, 1)
            k = minloc(fout(i:size(fout,1), 1), dim=1) + i - 1
            buf(:) = fout(i, :)
            fout(i, :) = fout(k, :)
            fout(k, :) = buf(:)
        end do
    end function core_basis
    
    function core_hamil(BIN) result(fout)
        integer, intent(in) :: BIN(:,:)
        real(dp), allocatable :: fout(:,:)
        integer i,j,a,sz
        integer nEd, nU, nDq, nup,ndn, nJ, nEp, nEc
        integer is(5), fs(5), ks(5), d(5)
        
        sz = size(BIN, 1)
        allocate(fout(sz,sz))
        
        fout(:,:) = 0.0
        
        do i = 1, sz
            do j = 1, sz
                if (i .eq. j) then
                    d = BIN(i, :)
                    nEd = d(1)
                    nU  = nEd*(nEd-1)/2
                    nDq = 6*popcnt(iand(d(2), 99)) - 4*popcnt(iand(d(2), 924))
                    nup = popcnt(iand(d(2), 992))
                    ndn = popcnt(iand(d(2), 31))
                    nJ  = (nup*(nup-1))/2 + (ndn*(ndn-1))/2
                    nEp = popcnt(d(3))
                    nEc = popcnt(d(4)) - popcnt(d(5))
                    
                    fout(i,i) = nEd*PMED - nEp*PMEP - nEc*PMEC + &
                                nU*PMU + nDq*PMDQ - nJ*PMJ
                else
                    is = BIN(i, :)
                    fs = BIN(j, :)
                    if ((fs(1) - is(1)) .eq. 1) then
                        a = ieor(is(2), fs(2))
                        if (popcnt(a) .eq. 1) then
                            ks = fs - is
                            if ((ks(3) .eq. a) .and. (ks(4) .eq. 0) .and. (ks(5) .eq. 0)) then
                                fout(i,j) = PMTLIST(wbit(a))
                            end if
                            if ((ks(3) .eq. 0) .and. (ks(4) .eq. a) .and. (ks(5) .eq. 0)) then
                                fout(i,j) = PMTC
                            end if
                            if ((ks(3) .eq. 0) .and. (ks(4) .eq. 0) .and. (ks(5) .eq. -a)) then
                                fout(i,j) = PMTC
                            end if
                        end if
                    end if
                    fout(j,i) = fout(i,j)
                end if
            end do
        end do
    end function core_hamil
    
end module modcore