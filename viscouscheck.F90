! $HeadURL$
! $Id$
PROGRAM viscouscheck
! Checks the error for the centered finite difference routines
! against the field sin(x) with known second derivative -sin(x).

 USE doublePrecision
 IMPLICIT NONE
 INTEGER :: n
 REAL(KIND = dp), PARAMETER :: pi = 4._dp*ATAN(1._dp)
 REAL(KIND = dp) :: h, hi
 REAL(KIND = dp) :: maxerr, err
 REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: x, u, d2u

 CALL get_command_argument (1, str)
 READ (str, fmt = '(I10)') n
 ALLOCATE ( x(0:n), u(0:n), d2u(0:n) )

! Setup grid x \in [0,2pi]
 h  = 2_dp*pi / n
 hi = 1_dp / h
 DO i = 0, n, 1
  x(i) = REAL(i,dp) * h
 END DO

 CALL printdble (x, n, 'x.txt')
 DO i = 0, n, 1
  u(i) = SIN(2_dp * pi * x(i))
 END DO

 up = 0_dp
 CALL 
 ! FIXME STARTHERE

#if WENOORDER == 3
#define RECONSTRUCT_FUNCTION reconstruct3
#elif WENOORDER == 5
#define RECONSTRUCT_FUNCTION reconstruct5
#else
  #error "WENOORDER not #defined or unknown"
#endif
  PRINT '(" WENO reconstruction order = ", I2)', WENOORDER

#ifdef VISCOUSORDER
  PRINT '(" Viscous term order = ", I2, " with viscosity = ", F16.12)', &
        VISCOUSORDER, nu
#if VISCOUSORDER == 2
#define VISCOUS_FUNCTION viscous2
#elif VISCOUSORDER == 4
#define VISCOUS_FUNCTION viscous4
#else
  #error "VISCOUSORDER unknown"
#endif
#else
#define VISCOUS_FUNCTION viscousnop
  PRINT '(" No viscous term present ")'
#endif

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
  CALL RECONSTRUCT_FUNCTION (f, u, n, -1)
  CALL printdble (f, n, 'urm.txt')
  CALL RECONSTRUCT_FUNCTION (f, u, n, +1)
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
  CALL RECONSTRUCT_FUNCTION (frp, fp, n, 1)
  CALL RECONSTRUCT_FUNCTION (frm, fm, n, -1)
  f = frp + frm
  CALL rhside(fp, f, n, hi)
  up(:,1) = u + dt * fp
  CALL VISCOUS_FUNCTION (up(:,1), u, n, hi, dt*nu)
! Substep 2
  CALL flux (fp, fm, up(:,1), n)
  CALL RECONSTRUCT_FUNCTION (frp, fp, n, 1)
  CALL RECONSTRUCT_FUNCTION (frm, fm, n, -1)
  f = frp + frm
  CALL rhside (fp, f, n, hi)
  up(:,2) = rk(1,1) * u + rk(1,2) * (up(:,1) + dt * fp)
  CALL VISCOUS_FUNCTION (up(:,2), up(:,1), n, hi, rk(1,2)*dt*nu)
! Substep 3
  CALL flux (fp, fm, up(:,2), n)
  CALL RECONSTRUCT_FUNCTION (frp, fp, n, 1)
  CALL RECONSTRUCT_FUNCTION (frm, fm, n, -1)
  f = frp + frm
  CALL rhside(fp, f, n, hi)
  u = rk(2,1) * u + rk(2,2) * (up(:,2) + dt * fp)
  CALL VISCOUS_FUNCTION (u, up(:,2), n, hi, rk(2,2)*dt*nu)
! Print solution every tprint step
  IF ((nt == nsteps) .OR. (MODULO(nt, tprint) == 0)) THEN
   WRITE (itstr, "(I7.7)") nt
   CALL printdble (u, n, 'sol' // TRIM(itstr) // '.txt')
  END IF
 END DO

 DEALLOCATE (x, xh, u, up, frp, frm, f, fp, fm)
 CALL cpu_time (finish)
 PRINT '(" CPU Time = ", f10.3, " seconds")', finish-start
END PROGRAM viscouscheck
