! $HeadURL$
! $Id$
SUBROUTINE viscousnop (vio, u, n, hi, alpha)
! This routine is a NOP with an interface matching to viscous{2,4}.
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n
 REAL(KIND = dp), INTENT(IN)     :: u(0:n), hi
 REAL(KIND = dp), INTENT(INOUT) :: vio(0:n)    ! where to accumulate result
 REAL(KIND = dp), INTENT(IN)     :: alpha      ! accumulation coefficient

END SUBROUTINE viscousnop
