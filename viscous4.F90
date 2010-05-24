! $HeadURL$
! $Id$
SUBROUTINE viscous4 (vio, u, n, hi, alpha)
! This accumulates $alpha * \frac{\partial^2}{\partial{}x^2} u$
! into vio using a fourth-order, centered finite difference.
 USE doublePrecision
 IMPLICIT NONE
 INTEGER,         INTENT(IN)    :: n
 REAL(KIND = dp), INTENT(IN)    :: u(0:n), hi
 REAL(KIND = dp), INTENT(INOUT) :: vio(0:n)   ! where to accumulate result
 REAL(KIND = dp), INTENT(IN)    :: alpha      ! accumulation coefficient
 REAL(KIND = dp)                :: c(0:2)     ! finite diff coefficients
 INTEGER                        :: i

 ! Precompute coefficients for (i-2), (i-1), (i), (i+1), (i+2) locations
 ! These incorporate the grid size, finite different weights, and alpha
 c(0) = alpha * (-30_dp/12_dp)*hi**2
 c(1) = alpha * (+16_dp/12_dp)*hi**2
 c(2) = alpha * (- 1_dp/12_dp)*hi**2

! Accumulate result in vio in almost a single pass
 vio(0) = vio(0) &
        + c(2)*u(n-2) + c(1)*u(n-1) + c(0)*u(0) + c(1)*u(1) + c(2)*u(2)
 vio(1) = vio(1) &
        + c(2)*u(n-1) + c(1)*u(0)   + c(0)*u(1) + c(1)*u(2) + c(2)*u(3)
 DO i = 2, n-2, 1
  vio(i) = vio(i) &
         + c(2)*u(i-2) + c(1)*u(i-1) + c(0)*u(i) + c(1)*u(i+1) + c(2)*u(i+2)
 END DO
 vio(n-1) = vio(n-1) &
          + c(2)*u(n-3) + c(1)*u(n-2) + c(0)*u(n-1) + c(1)*u(n) + c(2)*u(1)
 vio(n)   = vio(n) &
          + c(2)*u(n-2) + c(1)*u(n-1) + c(0)*u(n)   + c(1)*u(1) + c(2)*u(2)

END SUBROUTINE viscous4
