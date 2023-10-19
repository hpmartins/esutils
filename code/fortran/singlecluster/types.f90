module types
    implicit none

    integer, parameter :: sp = kind( 1.0 )
    integer, parameter :: dp = kind( 1.0d0 )

    integer, parameter :: wp = sp

    real(dp), parameter :: pi    = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: e_    = 2.7182818284590452353602874713527_dp
    complex(dp), parameter :: i_ = (0, 1)
    
    type imat
        integer, allocatable :: v(:,:)
    end type imat

    type rmat
        real(dp), allocatable :: v(:,:)
    end type rmat

    type dmat
        real(dp), allocatable :: v(:,:)
    end type dmat

    type ivec
        integer, allocatable :: v(:)
    end type ivec

    type rvec
        real(dp), allocatable :: v(:)
    end type rvec

    type dvec
        real(dp), allocatable :: v(:)
    end type dvec
end module
