! $HeadURL$
! $Id$
PROGRAM viscouscheck
! Checks the error for the centered finite difference routines
! against the field sin(x) with known second derivative -sin(x).

 USE doublePrecision
 IMPLICIT NONE
 REAL(KIND = dp), PARAMETER                   :: pi = 4._dp*ATAN(1._dp)
 INTEGER, PARAMETER                           :: ncases = 3
 INTEGER                                      :: powfinal, p, n, j
 REAL(KIND = dp)                              :: h
 REAL(KIND = dp), DIMENSION(:), ALLOCATABLE   :: u, expected
 REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: d2u, d2u_error
 CHARACTER (len = 20) :: str

 CALL get_command_argument (1, str)
 READ (str, fmt = '(I10)') powfinal
 ALLOCATE ( u(0:2**powfinal), expected(0:2**powfinal) )
 ALLOCATE ( d2u(ncases,0:2**powfinal), d2u_error(ncases,0:2**powfinal) )

 DO p = 3, powfinal, 1
   ! Initial manufactured field and clear accumulation locations
   n = 2**p
   CALL problem1(n, h, u, expected)
!  CALL printdble(u,        n, 'u.txt')         ! DEBUG
!  CALL printdble(expected, n, 'expected.txt')  ! DEBUG

   ! Compute the second derivatives of the field
   d2u = 0
   CALL viscousnop(d2u(1,:), u, n, 1_dp/h, REAL(1,dp))
   CALL viscous2(  d2u(2,:), u, n, 1_dp/h, REAL(1,dp))
   CALL viscous4(  d2u(3,:), u, n, 1_dp/h, REAL(1,dp))

!  CALL printdble(d2u(1,:),  n, 'd2u_1.txt')    ! DEBUG
!  CALL printdble(d2u(2,:),  n, 'd2u_2.txt')    ! DEBUG
!  CALL printdble(d2u(3,:),  n, 'd2u_3.txt')    ! DEBUG

   ! Compute pointwise errors
   d2u_error = d2u
   DO j = 1, ncases, 1
     d2u_error(j,:) = ABS(d2u_error(j,:) - expected)
   END DO

   ! Output results
   WRITE (*, "('n=',i8)", advance="no") n
   DO j = 1, ncases, 1
     WRITE (*, "(' ', g24.16)", advance="no") MAXVAL(d2u_error(j,:))
   END DO
   WRITE (*, "('')")

 END DO

 DEALLOCATE ( u, d2u, d2u_error, expected )
END PROGRAM viscouscheck

! For 1D grid size n, load sample data with spacing h into u and d2u
SUBROUTINE problem1 (n, h, u, d2u)
 USE doublePrecision
 IMPLICIT NONE
 REAL(KIND = dp), PARAMETER   :: pi = 4._dp*ATAN(1._dp)
 INTEGER, INTENT(IN)          :: n
 REAL(KIND = dp), INTENT(OUT) :: h, u(0:n), d2u(0:n)
 INTEGER                      :: i
 REAL(KIND = dp)              :: x

 h = 2_dp*pi/n
 DO i = 0, n, 1
  x      =  REAL(i,dp)*h
  u(i)   =  SIN(x)
  d2u(i) = -SIN(x)
 END DO
END SUBROUTINE problem1

! For 1D grid size n, load sample data with spacing h into u and d2u
SUBROUTINE problem2 (n, h, u, d2u)
 USE doublePrecision
 IMPLICIT NONE
 REAL(KIND = dp), PARAMETER   :: pi = 4._dp*ATAN(1._dp)
 INTEGER, INTENT(IN)          :: n
 REAL(KIND = dp), INTENT(OUT) :: h, u(0:n), d2u(0:n)
 INTEGER                      :: i
 REAL(KIND = dp)              :: x

 h = 2_dp*pi/n
 DO i = 0, n, 1
  x      =  REAL(i,dp)*h
  u(i)   =  COS(x + 2_dp*COS(3_dp*x))
  d2u(i) =   -(6_dp*SIN(3_dp*x) - 1_dp)**2 * COS(x + 2_dp*COS(3_dp*x)) &
            + 18_dp*SIN(x + 2_dp*COS(3_dp*x))*COS(3_dp*x)
 END DO
END SUBROUTINE problem2
