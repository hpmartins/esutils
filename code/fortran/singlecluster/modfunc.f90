module modfunc
use types
use moduniinv, only: uniinv
implicit none

contains

    subroutine basissiz(siz, st, coh)
        integer, intent(in) :: st, coh
        integer :: siz
        integer :: nst
        nst = popcnt(st)
        siz = 2**(10 - nst) + popcnt(iand(st, coh)) + &
              popcnt(iand(not(st), coh))*(2**(10 - nst - 1))
    end subroutine basissiz

    function vec2bin(a) result(fout) ! transforma array em binario
        integer, intent(in) :: a(10)
        integer :: fout
        integer :: c
        fout = 0
        do c = 1, 10
            fout = fout + a(c)*(2**(10 - c))
        end do
    end function vec2bin

    function wbits(a) result(fout)
        integer i, j
        integer, intent(in) :: a
        integer, allocatable :: fout(:)
        allocate(fout(popcnt(a)))
        j = 1
        do i = 0, 9
            if (btest(a, i)) then
                fout(j) = i
                j = j + 1
            end if
        end do
    end function wbits
    
    function wbit(a) result(fout)
        integer, intent(in) :: a
        integer :: i, fout
        
        fout = -1
        do i = 0, 9
            if (btest(a, i)) fout = i
        end do
    end function wbit 
    
    function tolower(str) result(fout)
        character(*), intent(in) :: str
        character(len(str)) :: fout
        integer :: i
        
        fout = str
        
        do i = 1, len(str)
            select case(str(i:i))
                case("A":"Z")
                    fout(i:i) = achar(iachar(str(i:i))+32)
            end select
        end do
    end function tolower
    
    function cfgname(b) result(fout)
        integer, intent(in) :: b(5)
        character(10) :: fout
        character(2) :: d, c, L1, L2, C1, C2
        
        write(d, '(i0)') b(1)

        if (popcnt(b(3)) .gt. 0) then
            L1 = 'L'
            if (popcnt(b(3)) .gt. 1) then
                write(L2, '(i0)') popcnt(b(3))
            else
                L2 = ''
            end if
        else
            L1 = ''
        end if
        
        if (popcnt(b(4)) .gt. 0) then
            C1 = 'C'
            if (popcnt(b(4)) .gt. 1) then
                write(C2, '(i0)') popcnt(b(4))
                C2 = trim(C2)
            else
                C2 = ''
            end if
        else
            C1 = ''
        end if
        
        if (popcnt(b(5)) .gt. 0) then
            c =  'c'
        else
            c = ''
        end if
        
        fout = trim('d'//trim(d)//trim(L1)//trim(L2)//trim(C1)//trim(C2)//trim(c))
    end function
    
    function occupation(b, v) result(fout)
        integer, intent(in) :: b(:,:)
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: bn(:,:), vn(:,:)
        real(dp) :: fout, tmp(1,1)
        allocate(bn(size(b, 1), 1))
        allocate(vn(size(v, 1), 1))
        bn(:, 1) = b(:, 1)
        vn(:, 1) = v(:)
        tmp = matmul(transpose(vn), bn*vn)
        fout = tmp(1,1)
    end function
    
    function i_grpsum(A, col) result(fout)
        integer, intent(in) :: A(:,:), col
        integer, allocatable :: fout(:,:), x(:), idx(:)
        integer :: i, j
        
        allocate(x(size(A, 1)))
        allocate(idx(size(A, 1)))
        
        x = A(:, col)
        call uniinv(x, idx)
        
        allocate(fout(maxval(idx), size(A, 2)))
        fout(:,:) = 0
        
        do j = 1, size(A, 2)
            do i = 1, maxval(idx)
                if (j .eq. col) then
                    fout(i,j) = minval(A(:,j), 1, idx .eq. i)
                else
                    fout(i,j) = sum(A(:,j), 1, idx .eq. i)
                end if
            end do
        end do
    end function i_grpsum
    
    function d_grpsum(A, col) result(fout)
        real(dp), intent(in)     :: A(:,:)
        integer,  intent(in)     :: col
        real(dp), allocatable    :: fout(:,:), x(:)
        integer,  allocatable    :: idx(:)
        integer :: i, j
        
        allocate(x(size(A, 1)))
        allocate(idx(size(A, 1)))
        
        x = A(:, col)
        call uniinv(x, idx)
        
        allocate(fout(maxval(idx), size(A, 2)))
        fout(:,:) = 0
        
        do j = 1, size(A, 2)
            do i = 1, maxval(idx)
                if (j .eq. col) then
                    fout(i,j) = minval(A(:,j), 1, idx .eq. i)
                else
                    fout(i,j) = sum(A(:,j), 1, idx .eq. i)
                end if
            end do
        end do
    end function d_grpsum
    
    function filtertol(A, col, tol) result(fout)
        real(dp), intent(in)     :: A(:,:), tol
        integer,  intent(in)     :: col
        real(dp), allocatable    :: fout(:,:)
        integer                  :: i
        allocate(fout(count(A(:,col) .gt. tol), size(A, 2)))
        if (size(A, 1) .ne. 0) then
            do i = 1, size(A, 2)
                fout(:, i) = pack(A(:, i), A(:,col) .gt. tol)
            end do
        else
            fout(:,:) = 0.0
        end if
    end function filtertol
    
    function findnl(s) result(fout)
        character(*), intent(in) :: s
        integer                  :: fout, i
        
        fout = len(s) + 1
        do i = 1, len(s)
            if (s(i:i) .eq. achar(10)) then
                fout = i
                return
            end if
        end do
    end function findnl
    
    function trim2(s) result(fout)
        character(len=*)          :: s
        character(:), allocatable :: fout
        integer                   :: i
        i = len_trim(s)
        fout = s(1:i)
    end function trim2
end module modfunc