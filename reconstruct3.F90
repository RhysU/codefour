! $HeadURL$
! $Id$
SUBROUTINE reconstruct3 (ur, u, n, bias)
! This function takes values of u in x0,x1,x2,x3,...,xn
! and returns 3rd-order WENO reconstructed values, ur, in x(1/2), x(3/2),...,x(n+1/2)
!
! This routine follows Liu, Osher, and Chan's 1994 paper section 3.5,
! and equation numbers refer to it.
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n, bias
 REAL(KIND = dp), INTENT(IN) :: u(0:n)
 REAL(KIND = dp), INTENT(OUT) :: ur(0:n)
 REAL(KIND = dp), PARAMETER :: eps = 1.d-6  ! guarantee nonzero denominator
 REAL(KIND = dp) :: beta(1:2)               ! smoothness indicators
 REAL(KIND = dp) :: w(1:2)                  ! nonlinear weights
 REAL(KIND = dp) :: wt(1:2), wtsumi         ! temporary nonlinear weights
 REAL(KIND = dp) :: urloc(1:2)              ! the two local reconstructions
 REAL(KIND = dp) :: a(1:2,1:2)              ! weights in reconstruction
 INTEGER :: i
 REAL(KIND = dp) :: v(-1:n+2)               ! add on periodic boundary conditions
 REAL(KIND = dp) :: v0, vp, vm              ! local values

 a(1,1) = -1.d0 / 2.d0
 a(1,2) =  3.d0 / 2.d0
 a(2,1) =  1.d0 / 2.d0
 a(2,2) =  1.d0 / 2.d0

! add on periodic boundary condition
! this is wasteful but results in a single loop so the code is easier to read
 v(0:n) = u(0:n)
 v(-1)  = u(n-1)
 v(n+1:n+2) = u(1:2)

 IF (bias > 0) THEN ! Bias to the left, case 1 in section 3.5
  DO i = 0, n, 1
   v0 = v(i)
   vp = v(i+1)
   vm = v(i-1)
! The reconstructed values at x(i+1/2) per p'(j), p'(j+1) from bottom of p205
! Note mistake in the p'j formula, i.e. (x-x).
   urloc(1) = a(1,1) * vm + a(1,2) * v0
   urloc(2) = a(2,1) * v0 + a(2,2) * vp
! Smoothness indicators from p206 just above equation 3.16
   beta(1) = (v0 - vm)**2
   beta(2) = (vp - v0)**2
! Compute nonlinear weights (3.17a)
   wt(1) = 0.5d0 / ((eps + beta(1))**2)
   wt(2) = 1.0d0 / ((eps + beta(2))**2)
   wtsumi = 1.d0 / (wt(1) + wt(2))
   w(1) = wt(1) * wtsumi
   w(2) = wt(2) * wtsumi
! Finally reconstruct, formula (3.16)
   ur(i) = w(1) * urloc(1) + w(2) * urloc(2)
  END DO
 ELSE ! biased to the right, case 2 in section 3.5
  DO i = 1, n+1, 1
   v0 = v(i)
   vp = v(i+1)
   vm = v(i-1)
! The reconstructed values at x(i-1/2) per p'(j), p'(j+1) from bottom of p205
! Note mistake in the p'j formula, i.e. (x-x).
   urloc(1) = a(2,1) * vm + a(2,2) * v0
   urloc(2) = a(1,2) * v0 + a(1,1) * vp
! Smoothness indicators from p206 just above equation 3.16
   beta(1) = (v0 - vm)**2
   beta(2) = (vp - v0)**2
! Compute nonlinear weights (3.17a)
   wt(1) = 1.0d0 / ((eps + beta(1))**2)
   wt(2) = 0.5d0 / ((eps + beta(2))**2)
   wtsumi = 1.d0 / (wt(1) + wt(2))
   w(1) = wt(1) * wtsumi
   w(2) = wt(2) * wtsumi
! Finally reconstruct, formula (3.16)
   ur(i-1) = w(1) * urloc(1) + w(2) * urloc(2)
  END DO
 END IF
END SUBROUTINE reconstruct3
