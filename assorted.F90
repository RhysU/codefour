! $HeadURL$
! $Id$

!> Print a double precision one dimensional array to a flat text file.
!!
!! @param u   Data to be output
!! @param n   Number of data points to write
!! @param str Filename to use
SUBROUTINE printdble (u, n, str)
  USE doublePrecision
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND = dp), INTENT(IN) :: u(0:n)
  CHARACTER(LEN = *), INTENT(IN) :: str
  INTEGER :: i

  OPEN  (2, FILE = TRIM(str), STATUS = 'UNKNOWN')
  WRITE (2, fmt = '(E24.16)') (u(i),i=0,n,1)
  CLOSE (2)
END SUBROUTINE printdble
