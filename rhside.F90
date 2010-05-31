! $HeadURL$
! $Id$

!> Evaluate the right hand side of the inviscid equation's time evolution
!! \f$
!! \partial_{t} \bar{u}_j = -\frac{1}{h}\left[
!!   f\left(u\left(x_{j+1/2},t\right)\right)
!!   -
!!   f\left(u\left(x_{j+1/2},t\right)\right)
!! \right]\f$
!! assuming periodic boundary conditions.  See section 2 of Liu,
!! Osher, and Chan's 1994 JCP paper for more details.
!!
!! @param fout The evaluated right hand side
!! @param fin  The input data \f$f\left(u\left(x_{j+1/2}\right)\right)\f$
!!             for \f$j\in\left\{0,\dots,n\right\}\f$.  Usually
!!             this will be an approximation found through reconstruction.
!! @param n    Grid size
!! @param hi   \f$\frac{1}{h}\f$
SUBROUTINE rhside (fout, fin, n, hi)
  USE doublePrecision
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND = dp), INTENT(IN) :: fin(0:n), hi
  REAL(KIND = dp), INTENT(OUT) :: fout(0:n)

  fout(1:n) = -hi * (fin(1:n) - fin(0:n-1))
  fout(0)   = fout(n)
END SUBROUTINE rhside
