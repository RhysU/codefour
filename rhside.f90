SUBROUTINE rhside (fout, fin, n, hi)
! This evaluates the right side assuming periodic boundary conditions
 USE doublePrecision
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: n
 REAL(KIND = dp), INTENT(IN) :: fin(0:n), hi
 REAL(KIND = dp), INTENT(OUT) :: fout(0:n)
 INTEGER :: i
 fout(1:n) = -hi * (fin(1:n) - fin(0:n-1))
 fout(0) = fout(n)
 
END SUBROUTINE rhside
