! $HeadURL$
! $Id$
SUBROUTINE printdble (u, n, str)
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n
 REAL(KIND = dp), INTENT(IN) :: u(0:n)
 CHARACTER (LEN = *), INTENT(IN) :: str
 INTEGER :: i
 OPEN (2, FILE = TRIM(str), STATUS = 'UNKNOWN')
 DO i = 0, n, 1
  WRITE (2, fmt = '(E24.16)') u(i)
 END DO
 CLOSE (2)
END SUBROUTINE printdble
