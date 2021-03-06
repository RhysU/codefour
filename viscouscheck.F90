! $HeadURL$
! $Id$

!> Checks the error for the centered finite difference routines
!! against a field with a known second derivative.
PROGRAM viscouscheck

 USE doublePrecision
 IMPLICIT NONE
 REAL(KIND = dp), PARAMETER                   :: pi = 4._dp*ATAN(1._dp)
 INTEGER, PARAMETER                           :: ncases = 3
 INTEGER                                      :: powfinal, p, n, j
 REAL(KIND = dp)                              :: h
 REAL(KIND = dp), DIMENSION(:), ALLOCATABLE   :: u, expected
 REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: d2u, d2u_error
 CHARACTER (len = 20) :: str

 IF (command_argument_count() /= 1) THEN
  CALL get_command_argument(0, str)
  print *, "Usage: ", trim(str), " 2**npoints"
  CALL EXIT(1)
 END IF

 CALL get_command_argument (1, str)
 READ (str, fmt = '(I10)') powfinal
 ALLOCATE ( u(0:2**powfinal), expected(0:2**powfinal) )
 ALLOCATE ( d2u(ncases,0:2**powfinal), d2u_error(ncases,0:2**powfinal) )

 DO p = 3, powfinal, 1
   ! Initial manufactured field and clear accumulation locations
   n = 2**p
   CALL problem2(n, h, u, expected)

   ! Compute the second derivatives of the field
   d2u = 0
   CALL viscousnop(d2u(1,:), u, n, 1_dp/h, REAL(1,dp))
   CALL viscous2(  d2u(2,:), u, n, 1_dp/h, REAL(1,dp))
   CALL viscous4(  d2u(3,:), u, n, 1_dp/h, REAL(1,dp))

   ! Compute pointwise errors
   d2u_error = d2u
   DO j = 1, ncases, 1
     d2u_error(j,:) = ABS(d2u_error(j,:) - expected)
   END DO

   ! Output results
   WRITE (*, "(i8)", advance="no") n
   DO j = 1, ncases, 1
     WRITE (*, "(' ', g24.16)", advance="no") MAXVAL(d2u_error(j,:))
   END DO
   WRITE (*, "('')")

 END DO

 DEALLOCATE ( u, d2u, d2u_error, expected )
END PROGRAM viscouscheck

!> Load data for the 1D sample problem \f$u(x) = \sin(x)\f$.
!!
!! @param n   Grid size
!! @param h   Uniform grid spacing
!! @param u   \f$u(x)\f$
!! @param d2u \f$\frac{\partial^2}{\partial{}x^2}u(x)\f$
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

!> Load data for the 1D sample problem
!! \f$u(x) = \cos\left(x+2\cos(3x)\right)\f$.
!!
!! @param n   Grid size
!! @param h   Uniform grid spacing
!! @param u   \f$u(x)\f$
!! @param d2u \f$\frac{\partial^2}{\partial{}x^2}u(x)\f$
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
  u(i)   =  COS(x + 2._dp*COS(3._dp*x))
  d2u(i) =   -(6._dp*SIN(3._dp*x) - 1._dp)**2 * COS(x + 2._dp*COS(3._dp*x)) &
            + 18._dp*SIN(x + 2._dp*COS(3._dp*x))*COS(3._dp*x)
 END DO
END SUBROUTINE problem2
