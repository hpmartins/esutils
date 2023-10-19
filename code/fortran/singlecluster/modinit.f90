module modinit
use types
implicit none

    ! Basis
    integer, allocatable :: BGND(:,:)
    integer, allocatable :: remlst(:), addlst(:)
    type(imat)           :: BREM(0:9), BADD(0:9)
    integer              :: gndsiz, remsizup, remsizdn, addsizup, addsizdn
    
    ! Hamiltonians
    real(dp), allocatable :: HGND(:,:), HCOR(:,:)
    type(rmat) :: HREM(0:9), HADD(0:9)
    
    ! Eigenvalues and eigenvectors
    real(dp) :: EVAGND
    real(dp), allocatable :: EVCGND(:), EVACOR(:), EVCCOR(:,:)
    type(dvec) :: EVAREM(0:9), EVAADD(0:9)
    type(dmat) :: EVCREM(0:9), EVCADD(0:9)
    
    ! Operators
    type(dmat) :: OPADDD(0:9), OPADDP(0:9), OPADDC(0:9)
    type(dmat) :: OPREMD(0:9), OPREMP(0:9), OPREMC(0:9)
    
    ! Transitions
    real(dp), allocatable :: ACOR(:,:)
    real(dp), allocatable :: AREMUP(:,:), AREMDN(:,:), AREMUPT(:,:), AREMDNT(:,:)
    real(dp), allocatable :: AADDUP(:,:), AADDDN(:,:), AADDUPT(:,:), AADDDNT(:,:)
    
    logical bool
    integer i,j,k,a,b,c,siz,st,narg
    character(256) EXCNM
    character(256) tmpstring

    !------------------------!
    !    input variables     !
    !------------------------!
    character(:), allocatable :: PMINPFL
    character(:), allocatable :: PMINPNM
    integer PMSTART ! start configuration
    real(dp) PMED, PMEP, PMDELTA, PMU, PMTS, PMTP, PMDQ, PMJ
    integer PMCOH   ! coherent position
    real(dp) PMDELTAC, PMTC, PMEC
    real(dp) PMTLIST(0:9)

end module