! $HeadURL$
! $Id$
SUBROUTINE reconstruct (ur, u, n, bias)
! This function takes values of u in x0,x1,x2,x3,...,xn
! and returns fifth-order WENO reconstructed values, ur, in x(1/2), x(3/2),...,x(n+1/2)
! 
! This routine follows Shu's 2009 SIAM Review paper and the equation numbers refer to it.
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n, bias
 REAL(KIND = dp), INTENT(IN) :: u(0:n)
 REAL(KIND = dp), INTENT(OUT) :: ur(0:n)
 REAL(KIND = dp), PARAMETER :: eps = 1.d-6  ! guarantee nonzero denominator
 REAL(KIND = dp) :: beta(1:3)               ! smoothness indicators
 REAL(KIND = dp) :: w(1:3)                  ! nonlinear weights
 REAL(KIND = dp) :: wt(1:3), wtsumi         ! temporary nonlinear weights
 REAL(KIND = dp) :: gam(1:3)                ! linear weights
 REAL(KIND = dp) :: urloc(1:3)              ! the three local reconstructions
 REAL(KIND = dp) :: a(1:3,1:3)              ! weights in reconstruction
 REAL(KIND = dp) :: b(1:2)                  ! constants for beta computation
 INTEGER :: i
 REAL(KIND = dp) :: v(-2:n+3)               ! add on periodic boundary conditions
 REAL(KIND = dp) :: v0, vp, vpp, vm, vmm    ! local values

 a(1,1) = 1.d0 / 3.d0
 a(1,2) = -7.d0 / 6.d0
 a(1,3) = 11.d0 / 6.d0
 a(2,1) = -1.d0 / 6.d0
 a(2,2) = 5.d0 / 6.d0
 a(2,3) = 1.d0 / 3.d0
 a(3,1) = 1.d0 / 3.d0
 a(3,2) = 5.d0 / 6.d0
 a(3,3) = -1.d0 / 6.d0

 b(1) = 13.d0 / 12.d0
 b(2) = 1.d0 / 4.d0
! just below (2.15)
 gam(1) = 1.d0 / 10.d0
 gam(2) = 3.d0 / 5.d0
 gam(3) = 3.d0 / 10.d0

! add on periodic boundary condition
! this is wasteful but results in a single loop so the code is easier to read
 v(0:n) = u(0:n)
 v(-2:-1) = u(n-2:n-1)
 v(n+1:n+3) = u(1:3)

 IF (bias > 0) THEN ! Bias to the left
  DO i = 0, n, 1
   v0 = v(i)
   vp = v(i+1)
   vpp = v(i+2)
   vm = v(i-1)
   vmm = v(i-2)
! The three reconstructed values at x(i+1/2)
! Formulas (2.11), (2.12), (2.13)
   urloc(1) = a(1,1) * vmm + a(1,2) * vm + a(1,3) * v0
   urloc(2) = a(2,1) * vm + a(2,2) * v0 + a(2,3) * vp
   urloc(3) = a(3,1) * v0 + a(3,2) * vp + a(3,3) * vpp
! Smoothness indicators, formula (2.17)
   beta(1) = b(1) * (vmm - 2.d0 * vm + v0)**2 + b(2) * (vmm - 4.d0 * vm + 3.d0 * v0)**2
   beta(2) = b(1) * (vm - 2.d0 * v0 + vp)**2 + b(2) * (vm - vp)**2
   beta(3) = b(1) * (v0 - 2.d0 * vp + vpp)**2 + b(2) * (3.d0 * v0 - 4.d0 * vp + vpp)**2
! Compute nonlinear weights (2.10)
   wt(1) = gam(1) / ((eps + beta(1))**2)
   wt(2) = gam(2) / ((eps + beta(2))**2)
   wt(3) = gam(3) / ((eps + beta(3))**2)
   wtsumi = 1.d0 / (wt(1) + wt(2) + wt(3))
   w(1) = wt(1) * wtsumi
   w(2) = wt(2) * wtsumi
   w(3) = wt(3) * wtsumi
! Finally reconstruct, formula (2.16)
   ur(i) = w(1) * urloc(1) + w(2) * urloc(2) + w(3) * urloc(3)
  END DO
 ELSE ! biased to the right
  DO i = 1, n+1, 1
   v0 = v(i)
   vp = v(i+1)
   vpp = v(i+2)
   vm = v(i-1)
   vmm = v(i-2)
! The three reconstructed values at x(i-1/2)
! Slightly different formulas than (2.11), (2.12), (2.13)
   urloc(1) = a(2,1) * vmm + a(2,2) * vm + a(2,3) * v0
   urloc(2) = a(3,1) * vm + a(3,2) * v0 + a(3,3) * vp
   urloc(3) = a(1,3) * v0 + a(1,2) * vp + a(1,1) * vpp
! Smoothness indicators, formula (2.17)
   beta(1) = b(1) * (vmm - 2.d0 * vm + v0)**2 + b(2) *( vmm - 4.d0 * vm + 3.d0 * v0)**2
   beta(2) = b(1) * (vm - 2.d0 * v0 + vp)**2 + b(2) * (vm - vp)**2
   beta(3) = b(1) * (v0 - 2.d0 * vp + vpp)**2 + b(2) * (3.d0 * v0 - 4.d0 * vp + vpp)**2
! Compute nonlinear weights (2.10)
   wt(1) = gam(3) / ((eps + beta(1))**2)
   wt(2) = gam(2) / ((eps + beta(2))**2)
   wt(3) = gam(1) / ((eps + beta(3))**2)
   wtsumi = 1.d0 / (wt(1) + wt(2) + wt(3))
   w(1) = wt(1) * wtsumi
   w(2) = wt(2) * wtsumi
   w(3) = wt(3) * wtsumi
! Finally reconstruct! Formula (2.16)
   ur(i-1) = w(1) * urloc(1) + w(2) * urloc(2) + w(3) * urloc(3)
  END DO
 END IF
END SUBROUTINE reconstruct
