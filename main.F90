! $HeadURL$
! $Id$
PROGRAM main
! This program solves the one-dimensional scalar conservation law
! u_t + f(u)_x =  0 using WENO reconstruction per reconstructWENO_ORDER 
! #define and a Lax-Friedrichs flux (in flux.f90).  Time stepping is performed
! using a 3rd-order TVD Runge-Kutta scheme.  Periodic boundary conditions are
! used, but that is hidden in subroutines.

  USE HDF5
  USE H5LT
  USE doublePrecision
  IMPLICIT NONE

! Problem parameters
  REAL(KIND = dp), PARAMETER :: pi      = 4._dp*ATAN(1._dp)
  REAL(KIND = dp), PARAMETER :: tend    = 0.30_dp
  REAL(KIND = dp), PARAMETER :: cfl     = 0.5_dp
  REAL(KIND = dp), PARAMETER :: nu      = 1._dp/500._dp
  REAL(KIND = dp), PARAMETER :: rk(2,2) = RESHAPE( &              ! TVD RK3
             [0.75_dp, (1._dp/3._dp), (0.25_dp), (2._dp/3._dp)], SHAPE(rk))

! Working storage
  INTEGER :: i, n, nt, nsteps
  REAL(KIND = dp) :: h, hi, t, dt
  REAL(KIND = dp), DIMENSION(:),   ALLOCATABLE :: x, xh, u, frp, frm, fp, fm, f
  REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: up
  CHARACTER(LEN = 20) :: str

! Handles related to to HDF5 operations
  INTEGER(HID_T) :: file_hid, filespace_hid, dataspace_hid, dataset_hid
  INTEGER(HID_T) :: memspace_hid
  INTEGER        :: error

! Sanity check incoming argument count
  IF (command_argument_count() /= 2) THEN
    CALL get_command_argument(0, str)
    print *, "Usage: ", trim(str), " outputfile npoints"
    CALL EXIT(1)
  END IF

! Open the output file from arguments
  CALL h5open_f(error)
  CALL get_command_argument (1, str)
  CALL h5fcreate_f(trim(str), H5F_ACC_TRUNC_F, file_hid, error)

! Determine the grid size from arguments and allocate storage
  CALL get_command_argument (2, str)
  READ (str, fmt = '(I10)') n
  ALLOCATE ( x(0:n), xh(0:n), u(0:n), up(0:n,1:2), &
             frp(0:n), frm(0:n), f(0:n), fp(0:n), fm(0:n) )

! Setup grid x \in [0,1] containing n points
  h = 1.d0 / n
  hi = 1.d0 / h
  FORALL(i=0:n:1) x(i) = i*h
  xh = x + 0.5_dp*h
  CALL h5ltmake_dataset_double_f( &
       file_hid, "x", 1, [INTEGER(HSIZE_T)::n], x, error)
  CALL h5ltmake_dataset_double_f( &
       file_hid, "xh", 1, [INTEGER(HSIZE_T)::n], xh, error)

! Determine the WENO reconstruction order from #defines
#if WENOORDER == 3
#define RECONSTRUCT_FUNCTION reconstruct3
#elif WENOORDER == 5
#define RECONSTRUCT_FUNCTION reconstruct5
#else
  #error "WENOORDER not #defined or unknown"
#endif
  PRINT '(" WENO reconstruction order = ", I2)', WENOORDER
  CALL h5ltmake_dataset_int_f( &
       file_hid, "weno", 1, [INTEGER(HSIZE_T)::1], [WENOORDER], error)

! Determine viscous term presence and/or order from #defines
#ifdef VISCOUSORDER
  PRINT '(" Viscous term order = ", I2, " with viscosity = ", F16.12)', &
        VISCOUSORDER, nu
  CALL h5ltmake_dataset_int_f( &
       file_hid, "viscous", 1, [INTEGER(HSIZE_T)::1], [VISCOUSORDER], error)
  CALL h5ltmake_dataset_double_f( &
       file_hid, "nu", 1, [INTEGER(HSIZE_T)::1], [nu], error)
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

! Determine stable time step size
  dt = cfl * h  ! Convective stability
#ifdef VISCOUSORDER
! Constant from 1 + x + x**2/2 + x**3/6 ~ 1 along real axis per TVD RK3
  dt = MIN(dt, 2.512_dp / (nu * pi**2 * n **2)) ! Diffusive stability
#endif
  nsteps = INT(tend / dt)
  dt = tend / REAL(nsteps, dp)
  PRINT '(" Number of time steps = ", I7, " with dt = ", F16.12)', &
        nsteps, dt
  CALL h5ltmake_dataset_double_f( &
       file_hid, "t", 1, [INTEGER(HSIZE_T)::nsteps+1], &
       [(dt*i, i=0, nsteps, 1)], error)

! Initial condition 1: simple sine wave for debugging
! u = SIN(2._dp*pi*x)

! Initial condition 2: viscous analytic solution from Hopf-Cole
  t = 0
  u =   (                        27*pi*SIN(1 + 6*pi*x) ) &
      / ( 125*(10*EXP(9*pi*pi*t/25)+ 9*COS(1 + 6*pi*x)))

! Create dataset to store initial condition and (space x time) solution
  CALL h5screate_simple_f( &
       2, [INTEGER(HSIZE_T)::n, nsteps + 1], dataspace_hid, error)
  CALL h5dcreate_f( &
       file_hid, "u", H5T_NATIVE_DOUBLE, dataspace_hid, dataset_hid, error)
! Create dataspace describing how solution is stored in memory
  CALL h5screate_simple_f(2, [INTEGER(HSIZE_T)::n, 1], memspace_hid, error)
! Create dataspace describing how one solution time is stored on disk
  CALL h5dget_space_f(dataset_hid, filespace_hid, error)

! Write the initial condition
  CALL h5sselect_hyperslab_f(filespace_hid, H5S_SELECT_SET_F, &
       [INTEGER(HSIZE_T)::0, 0], [INTEGER(HSIZE_T)::n, 1], error)
  CALL h5dwrite_f(dataset_hid, H5T_NATIVE_DOUBLE, u(0:n-1), &
       [INTEGER(HSIZE_T)::n, 1], error, memspace_hid, filespace_hid)

! Advance the solution in time using explicit TVD RK3
  timeloop: DO nt = 1, nsteps, 1
    t = REAL(nt-1,dp) * dt

!   Substep 1
    CALL flux (fp, fm, u, n)
    CALL RECONSTRUCT_FUNCTION (frp, fp, n, 1)
    CALL RECONSTRUCT_FUNCTION (frm, fm, n, -1)
    f = frp + frm
    CALL rhside(fp, f, n, hi)
    CALL VISCOUS_FUNCTION (fp, u, n, hi, nu)
    up(:,1) = u + dt * fp

!   Substep 2
    CALL flux (fp, fm, up(:,1), n)
    CALL RECONSTRUCT_FUNCTION (frp, fp, n, 1)
    CALL RECONSTRUCT_FUNCTION (frm, fm, n, -1)
    f = frp + frm
    CALL rhside (fp, f, n, hi)
    CALL VISCOUS_FUNCTION (fp, up(:,1), n, hi, nu)
    up(:,2) = rk(1,1) * u + rk(1,2) * (up(:,1) + dt * fp)

!   Substep 3
    CALL flux (fp, fm, up(:,2), n)
    CALL RECONSTRUCT_FUNCTION (frp, fp, n, 1)
    CALL RECONSTRUCT_FUNCTION (frm, fm, n, -1)
    f = frp + frm
    CALL rhside(fp, f, n, hi)
    CALL VISCOUS_FUNCTION (fp, up(:,2), n, hi, nu)
    u = rk(2,1) * u + rk(2,2) * (up(:,2) + dt * fp)

!   Write the solution at the current time step
    CALL h5sselect_hyperslab_f(filespace_hid, H5S_SELECT_SET_F, &
         [INTEGER(HSIZE_T)::0, nt], [INTEGER(HSIZE_T)::n, 1], error)
    CALL h5dwrite_f(dataset_hid, H5T_NATIVE_DOUBLE, u(0:n-1), &
         [INTEGER(HSIZE_T)::n, 1], error, memspace_hid, filespace_hid)
  END DO timeloop

  ! Tear down resources
  DEALLOCATE (x, xh, u, up, frp, frm, f, fp, fm)
  CALL h5sclose_f(memspace_hid, error)
  CALL h5dclose_f(dataset_hid, error)
  CALL h5sclose_f(dataspace_hid, error)
  CALL h5fclose_f(file_hid, error)
  CALL h5close_f(error)

END PROGRAM main
