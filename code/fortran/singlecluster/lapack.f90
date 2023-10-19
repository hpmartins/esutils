module lapack
implicit none

! This is the precision that LAPACK "d" routines were compiled with (typically
! double precision, unless a special compiler option was used while compiling
! LAPACK). This "dp" is only used in lapack.f90
! The "d" routines data type is defined as "double precision", so
! we make "dp" the same kind as 0.d0 ("double precision"), so
! as long as LAPACK and this file were compiled with the same compiler options,
! it will be consistent. (If for example all double precision is promoted to
! quadruple precision, it will be promoted both in LAPACK and here.)
integer, parameter:: dp=kind(0.d0)

interface

    SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
                       LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * )
    END SUBROUTINE

end interface

contains

end module
