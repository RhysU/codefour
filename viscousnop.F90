! $HeadURL$
! $Id$

!> Do nothing. Provided for interface similarity to viscous2() and viscous4()
!! to support solving selectively solving inviscid equations using
!! compile-time definitions.
!!
!! @param v     Ignored
!! @param u     Ignored
!! @param n     Ignored
!! @param hi    Ignored
!! @param alpha Ignored
SUBROUTINE viscousnop (v, u, n, hi, alpha)
  USE doublePrecision
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND = dp), INTENT(IN)     :: u(0:n), hi
  REAL(KIND = dp), INTENT(INOUT) :: v(0:n)      ! where to accumulate result
  REAL(KIND = dp), INTENT(IN)     :: alpha      ! accumulation coefficient

END SUBROUTINE viscousnop
