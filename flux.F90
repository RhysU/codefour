! $HeadURL$
! $Id$

!> Compute the Lax-Friedrichs flux
!! \f$\hat{f}^{\text{LF}}\left(u^{-},u^{+}\right)\f$ given values of
!! \f$u\f$.  See section 3.1 of Shu's 2009 SIAM Review paper
!! or section 2 of Liu, Sher, and Chan's 1994 JCP paper for more details.
!!
!! @param fp
!! @param fm
!! @param u  \f$u\f$
!! @param n  Grid size
SUBROUTINE flux (fp, fm, u, n)
  USE doublePrecision
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: n
  REAL(KIND = dp), INTENT(IN)  :: u(0:n)
  REAL(KIND = dp), INTENT(OUT) :: fp(0:n), fm(0:n)
  REAL(KIND = dp)              :: alpha

! Euler equation where f(u) := u**2/2
  alpha = MAXVAL(ABS(u))
  fp    = 0.5_dp * (0.5_dp*u**2 + alpha*u)
  fm    = 0.5_dp * (0.5_dp*u**2 - alpha*u)

! The scalar advection equation where f(u) := u
! alpha = 1.d0
! fp    = 0.5_dp * (u + alpha*u)
! fm    = 0.5_dp * (u - alpha*u)

END SUBROUTINE flux
