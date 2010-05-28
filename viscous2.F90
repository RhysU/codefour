! $HeadURL$
! $Id$

!> Compute \f$v \leftarrow{} v + \alpha \frac{\partial^2}{\partial{}x^2} u\f$
!! using a second-order, centered finite difference.  The differencing requires
!! that the grid be uniformly spaced.  The differencing assumes the
!! grid is periodic with <tt>v(0) == v(n)</tt> and <tt>u(0) == u(n)</tt>.
!!
!! @param v     \f$v\f$
!! @param u     \f$u\f$
!! @param n     Grid size
!! @param hi    \f$1/h\f$ where \f$h\f$ is the uniform grid spacing
!! @param alpha Multiplicative coefficient \f$\alpha\f$
SUBROUTINE viscous2 (v, u, n, hi, alpha)
  USE doublePrecision
  IMPLICIT NONE
  INTEGER,         INTENT(IN)    :: n
  REAL(KIND = dp), INTENT(IN)    :: u(0:n), hi
  REAL(KIND = dp), INTENT(INOUT) :: v(0:n)     ! where to accumulate result
  REAL(KIND = dp), INTENT(IN)    :: alpha      ! accumulation coefficient
  REAL(KIND = dp)                :: c(0:1)     ! finite diff coefficients
  INTEGER                        :: i

! Precompute coefficients for (i-1), (i), (i+1) locations
! These incorporate the grid size, finite different weights, and alpha
  c(0) = alpha * (-2._dp)*hi**2
  c(1) = alpha * (+1._dp)*hi**2

! Accumulate result in v in almost a single pass
  v(0) = v(0) + c(1)*u(n-1) + c(0)*u(0) + c(1)*u(1)
  DO i = 1, n-1, 1
    v(i) = v(i) + c(1)*u(i-1) + c(0)*u(i) + c(1)*u(i+1)
  END DO
  v(n) = v(n) + c(1)*u(n-1) + c(0)*u(n) + c(1)*u(1)

END SUBROUTINE viscous2
