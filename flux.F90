! $HeadURL$
! $Id$
SUBROUTINE flux (fp, fm, u, n)
! This function takes values of u and computes the two fluxes
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n
 REAL(KIND = dp), INTENT(IN) :: u(0:n)
 REAL(KIND = dp), INTENT(OUT) :: fp(0:n), fm(0:n)
 REAL(KIND = dp) :: alpha
 INTEGER :: i

 IF (.FALSE.) THEN
! The scalar advection equation
  DO i = 0, n, 1
   alpha = 1.d0
   fp(i) = 0.5_dp * (u(i) + alpha * u(i))
   fm(i) = 0.5_dp * (u(i) - alpha * u(i))
  END DO
 END IF

 IF (.TRUE.) THEN
! Eulers' equation
  DO i = 0, n, 1
   alpha = MAXVAL(ABS(u))
   fp(i) = 0.5_dp * (0.5_dp * u(i)**2 + alpha * u(i))
   fm(i) = 0.5_dp * (0.5_dp * u(i)**2 - alpha * u(i))
  END DO
 END IF

END SUBROUTINE flux
