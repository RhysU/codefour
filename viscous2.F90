! $HeadURL$
! $Id$
SUBROUTINE viscous2 (vio, u, n, hi, alpha)
! This accumulates $alpha * \frac{\partial^2}{\partial{}x^2} u$
! into vio using a second-order, centered finite difference.
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n
 REAL(KIND = dp), INTENT(IN)  :: u(0:n), hi
 REAL(KIND = dp), INTENT(OUT) :: vio(0:n)      ! where to accumulate result
 REAL(KIND = dp), INTENT(IN)  :: alpha         ! accumulation coefficient
 REAL(KIND = dp)              :: c(0:1)        ! finite diff coefficients
 INTEGER                      :: i

! Precompute coefficients for (i-1), (i), (i+1) locations
! These incorporate the grid size, finite different weights, and alpha
 c(0) = alpha * (-2_dp)*hi**2
 c(1) = alpha * (+1_dp)*hi**2

! Accumulate result in vio in almost a single pass
 vio(0) = vio(0) + c(1)*u(n-1) + c(0)*u(0) + c(1)*u(1)
 DO i = 1, n-1, 1
  vio(i) = vio(i) + c(1)*u(i-1) + c(0)*u(i) + c(1)*u(i+1)
 END DO
 vio(n) = vio(n) + c(1)*u(n-1) + c(0)*u(n) + c(1)*u(1)

END SUBROUTINE viscous2
