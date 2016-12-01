 SUBROUTINE INIT_MESH_METRICS

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: IV

 WRITE(*, '(A)', ADVANCE = 'NO') ' INITIALIZING MESH METRICS... '

 !Allocate memory
 ALLOCATE(XVMP(NVERT), YVMP(NVERT))
 ALLOCATE(XVMP_VEL(NVERT), YVMP_VEL(NVERT))
 ALLOCATE(XXSolu(2,N,N,NCELL))
 ALLOCATE(Jac(N,N,NCELL), S1(2,3,N+1,N,NCELL), S2(2,3,N,N+1,NCELL))

 XVMP_VEL(:) = 0._DP; YVMP_VEL(:) = 0._DP

 DO IV = 1, NVERT
	XVMP(IV) = XV(IV)
	YVMP(IV) = YV(IV)
 END DO

 CALL CALCULATE_JACOBIAN

 CALL GetXXSolu

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''
 
 END SUBROUTINE INIT_MESH_METRICS
 !-------------------------------------------------------------------------------!


 SUBROUTINE UPDATE_MESH_METRICS

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER :: IV

 REAL(KIND=DP) :: tol
 REAL(KIND=DP) :: xd, yd, rad, theta
 REAL(KIND=DP) :: dr, vel_r
 REAL(KIND=DP) :: amp, x_L, y_L, t0

 IF (INIT_COND .EQ. SS_VORTEX) THEN

	amp = 0.05_DP; t0 = 1._DP

	DO IV = 1, NVERT
		xd = XV(IV); yd = YV(IV)
		rad   = SQRT(xd**2 + yd**2)

		IF ( (rad .GE. 1.05_DP) .AND. (rad .LE. 1.38_DP) .AND. &
			 (xd .GT. amp) .AND. (yd .GT. amp) ) THEN

			theta = ATAN(yd/xd)
			dr    = amp*SIN(4._DP*pi*theta/(0.5_DP*pi))*SIN(4._DP*pi*ctime/t0)
			vel_r = 4._DP*pi/t0*amp*SIN(4._DP*pi*theta/(0.5_DP*pi))*COS(4._DP*pi*ctime/t0)

			XVMP(IV) = (rad+dr)*COS(theta)
			YVMP(IV) = (rad+dr)*SIN(theta)

			XVMP_VEL(IV) = vel_r*COS(theta)
			YVMP_VEL(IV) = vel_r*SIN(theta)

		END IF
	END DO

 ELSEIF (INIT_COND .EQ. X_SINWAVE .OR. INIT_COND .EQ. FREE_STRM) THEN

	amp = 0.075_DP
	x_L = 1._DP; y_L = 1._DP; t0 = 1._DP
	tol = 1.0E-4

	DO IV = 1, NVERT
		IF ((XV(IV) .GT. -1.+tol) .AND. (XV(IV) .LT. 1.-tol) .AND. &
			(YV(IV) .GT. -1.+tol) .AND. (YV(IV) .LT. 1.-tol)) THEN

			xd = XV(IV); yd = YV(IV)

			XVMP(IV) = xd + amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*SIN(2._DP*pi*ctime/t0)
			YVMP(IV) = yd + amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*SIN(4._DP*pi*ctime/t0)

			XVMP_VEL(IV) = 2._DP*pi/t0*amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*COS(2._DP*pi*ctime/t0)
			YVMP_VEL(IV) = 4._DP*pi/t0*amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*COS(4._DP*pi*ctime/t0)
		END IF
	END DO

 ELSEIF (INIT_COND .EQ. EU_VORTEX) THEN

	amp = 0.2_DP
	x_L = 5._DP; y_L = 5._DP; t0 = 1.5_DP
	tol = 1.0E-3

	DO IV = 1, NVERT
		IF ((XV(IV) .GT. -5.+tol) .AND. (XV(IV) .LT. 5.-tol) .AND. &
			(YV(IV) .GT. -5.+tol) .AND. (YV(IV) .LT. 5.-tol)) THEN

			xd = XV(IV); yd = YV(IV)

			XVMP(IV) = xd + amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*SIN(2._DP*pi*ctime/t0)
			YVMP(IV) = yd + amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*SIN(4._DP*pi*ctime/t0)

			XVMP_VEL(IV) = 2._DP*pi/t0*amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*COS(2._DP*pi*ctime/t0)
			YVMP_VEL(IV) = 4._DP*pi/t0*amp*SIN(pi*xd/x_L)*SIN(pi*yd/y_L)*COS(4._DP*pi*ctime/t0)

		END IF
	END DO

 END IF

 CALL CALCULATE_JACOBIAN

 CALL GetXXSolu

 END SUBROUTINE UPDATE_MESH_METRICS
 !-------------------------------------------------------------------------------!
