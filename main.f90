! $HeadURL$
! $Id$
PROGRAM main
! This program solves the one-dimensional scalar conservation law u_t + f(u)_x = 0
! using fifth-order WENO reconstruction and a Lax-Friedrichs flux (in flux.f90).
! Time stepping is performed using a 3rd-order TVD Runge-Kutta scheme.  Periodic
! boundary conditions are used, but that is hidden in subroutines. 

 USE doublePrecision
 IMPLICIT NONE
 INTEGER :: n
 REAL(KIND = dp), PARAMETER :: tend = 0.30d0
 REAL(KIND = dp), PARAMETER :: cfl = 0.5d0
 REAL(KIND = dp), PARAMETER :: pi = 3.141592653589793d0
 INTEGER, PARAMETER :: tprint = 1    
 INTEGER :: i, j, nt, nsteps
 REAL(KIND = dp) :: h, hi, t, dt, lambda, pi
 REAL(KIND = dp) :: maxerr, err, etmp, rk(1:2,1:2)
 REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: x, xh, u, frp, frm, fp, fm, f
 REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: up
 CHARACTER (len = 20) :: str
 CHARACTER (len = 100) :: itstr
 REAL :: start, finish

 CALL get_command_argument (1, str)
 READ (str, fmt = '(I10)') n
 CALL cpu_time (start)
 ALLOCATE (x(0:n), xh(0:n), u(0:n), up(0:n,1:2), frp(0:n), frm(0:n), f(0:n), fp(0:n), fm(0:n))
 h = 1.d0 / n
 hi = 1.d0 / h
! Setup grid x \in [0,1]
 DO i = 0, n, 1
  x(i) = REAL(i,dp) * h
  xh(i) = x(i) + 0.5_dp * h
 END DO
 CALL printdble (x, n, 'x.txt')
 CALL printdble (xh, n, 'xh.txt')
! Initial data
 u = 0_dp
 DO i = 0, n, 1
  u(i) = SIN(2_dp * pi * x(i))
! if ((x(i) > 0.25) .AND. (x(i) < 0.75)) u(i) = 1_dp
 END DO

! A reconstruction example
! uncomment to check the order of accuracy of the reconstruction
 IF (.FALSE.) THEN
  DO i = 0, n, 1
! Need averages
   u(i) = -1.d0 / (h * 20.d0 * pi) * ( &
      COS(20.d0 * pi * (x(i) + 0.5d0 * h)) - COS(20.d0 * pi * (x(i) - 0.5d0 * h)))
   fp(i) = SIN(20.d0 * pi * xh(i))
  END DO
  CALL printdble (fp, n, 'usin.txt')
  CALL reconstruct (f, u, n, -1)
  CALL printdble (f, n, 'urm.txt')
  CALL reconstruct (f, u, n, +1)
  CALL printdble (f, n, 'urp.txt')
  maxerr = 0
  DO i = 0, n, 1
   err = ABS(f(i) - fp(i))
   maxerr = maxerr + h * err**2
  END DO
  maxerr = SQRT(maxerr)
  WRITE (*, *) h, maxerr
 END IF

 t = 0_dp
 dt = cfl * h
 lambda = dt * hi
 nsteps = INT(tend / dt)
 dt = tend / REAL(nsteps, dp)
 PRINT '(" Number of time steps = ", I7)', nsteps
! TVD Runge-Kutta third-order accurate
 rk(1,1) = 0.75_dp
 rk(1,2) = 0.25_dp
 rk(2,1) = REAL(1,dp) / REAL(3,dp)
 rk(2,2) = REAL(2,dp) / REAL(3,dp)

! time loop
 DO nt = 1, nsteps, 1
  t = REAL(nt-1,dp) * dt
! Substep 1
  CALL flux (fp, fm, u, n)
  CALL reconstruct (frp, fp, n, 1)
  CALL reconstruct (frm, fm, n, -1)
  f = frp + frm 
  CALL rhside(fp, f, n, hi)
  up(:,1) = u + dt * fp
! Substep 2
  CALL flux (fp, fm, up(:,1), n)
  CALL reconstruct (frp, fp, n, 1)
  CALL reconstruct (frm, fm, n, -1)
  f = frp + frm
  CALL rhside (fp, f, n, hi)
  up(:,2) = rk(1,1) * u + rk(1,2) * (up(:,1) + dt * fp)
! Substep 3
  CALL flux (fp, fm, up(:,2), n)
  CALL reconstruct (frp, fp, n, 1)
  CALL reconstruct (frm, fm, n, -1)
  f = frp + frm
  CALL rhside(fp, f, n, hi)
  u = rk(2,1) * u + rk(2,2) * (up(:,2) + dt * fp)
! Print solution every tprint step
  IF ((nt == nsteps) .OR. (MODULO(nt, tprint) == 0)) THEN
   WRITE (itstr, "(I7.7)") nt
   CALL printdble (u, n, 'sol' // TRIM(itstr) // '.txt')
  END IF
 END DO

 DEALLOCATE (x, xh, u, up, frp, frm, f, fp, fm)
 CALL cpu_time (finish)
 PRINT '(" CPU Time = ", f10.3, " seconds")', finish-start
END PROGRAM main
