module linalg
  use types, only: dp
  use lapack, only: dsyevd
  use utils, only: stop_error
  implicit none
  private
  public eigh


  ! eigenvalue/-vector problem for real symmetric/complex hermitian matrices:
  interface eigh
     module procedure deigh_simple
  end interface eigh

  interface assert_shape
     module procedure dassert_shape
  end interface assert_shape
  
contains

  subroutine deigh_simple(Am, lam, c)
    ! solves eigen value problem for all eigenvalues and eigenvectors
    ! Am must by symmetric
    ! Only the lower triangular part of Am is used.
    real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.
    integer :: n
    ! lapack variables
    integer :: lwork, liwork, info
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)

    ! solve
    n = size(Am,1)
    call assert_shape(Am, [n, n], "eigh", "Am")
    call assert_shape(c, [n, n], "eigh", "c")
    lwork = 1 + 6*n + 2*n**2
    liwork = 3 + 5*n
    allocate(work(lwork), iwork(liwork))
    c = Am
    call dsyevd('V','L',n,c,n,lam,work,lwork,iwork,liwork,info)
    if (info /= 0) then
       print *, "dsyevd returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "the algorithm failed to compute an eigenvalue while working"
          print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
          print *, "through", mod(info, n+1)
       end if
       call stop_error('eigh: dsyevd error')
    end if
  end subroutine deigh_simple

  subroutine dassert_shape(A, shap, routine, matname)
    ! make sure a given real matrix has a given shape
    real(dp), intent(in) :: A(:,:)
    integer, intent(in) :: shap(:)
    character(len=*) :: routine, matname

    if(any(shape(A) /= shap)) then
       print *, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print *, "Shape should be ", shap
       call stop_error("Aborting due to illegal matrix operation")
    end if
  end subroutine dassert_shape
end module linalg
