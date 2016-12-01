 SUBROUTINE ITERATIONS

 USE setup
 USE MESH2D
 IMPLICIT NONE

 CHARACTER (LEN = 100) :: ERR_MSG

 !Allocate Runge-Kutta memory
 IF (RK_stage .EQ. 4) THEN
	ALLOCATE(RK_k1(4,N,N,NCELL), &
			 RK_k2(4,N,N,NCELL), &
			 RK_k3(4,N,N,NCELL), &
			 RK_k4(4,N,N,NCELL))
 ELSEIF (RK_stage .EQ. 5) THEN
	ALLOCATE(RK_k1(4,N,N,NCELL), &
			 RK_k2(4,N,N,NCELL), &
			 RK_k3(4,N,N,NCELL), &
			 RK_k4(4,N,N,NCELL), &
			 RK_k5(4,N,N,NCELL))
 END IF

 DO WHILE (ITER .LT. ITERMAX)

	!Call appropriate time integration subroutine
	SELECT CASE(RK_stage)
	CASE(1)
		CALL STEP_EULER
	CASE(4)
		CALL STEP_RK4
	CASE(5)
		CALL STEP_SSPRK5
	CASE DEFAULT
		WRITE(ERR_MSG,'(A,I2)') ' INVALID TIME INTEGRATION SCHEME ID: ', OUT_FORMAT
		CALL SHUTDOWN (ERR_MSG)
	END SELECT

	ITER = ITER + 1

	CALL check_min_density

	IF ((ITER .EQ. 1) .OR. (MOD(ITER,10) .EQ. 0) .OR. (ITER .EQ. ITERMAX)) THEN
		CALL resid_norm
		CALL comp_error
	END IF

	IF ((MOD(ITER,TECFREQ) .EQ. 0) .OR. (ITER .EQ. ITERMAX)) THEN
		SELECT CASE (OUT_FORMAT)
		CASE (T_ASCII)
			CALL WRITE_TEC_DATA_NODAL
		CASE (T_BINARY)
			CALL WRITE_TEC_DATA_NODAL_BINARY
		CASE DEFAULT
			WRITE(ERR_MSG,'(A,I2)') ' INVALID OUTPUT FORMAT: ', OUT_FORMAT
			CALL SHUTDOWN (ERR_MSG)
		END SELECT
	END IF

 END DO
  
 END SUBROUTINE ITERATIONS
!----------------------------------------------------------------------------!


 SUBROUTINE STEP_EULER

 USE setup
 USE MESH2D
 IMPLICIT NONE

 REAL(KIND=DP), DIMENSION(4,N,N,NCELL) :: Q0

 Q0 = Q

 CALL compresid
 Q = Q0 + dt*resid
 ctime = ctime + dt 
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 END SUBROUTINE STEP_EULER
 !----------------------------------------------------------------------------!


 SUBROUTINE STEP_RK4

 USE setup
 USE MESH2D
 IMPLICIT NONE

 REAL(KIND=DP), DIMENSION(4,N,N,NCELL) :: Q0
 REAL(KIND=DP) :: b1, b2, b3, b4
 REAL(KIND=DP) :: ctime_n

 !--------------------------------------------------------------------!
 !Define RK coefficients:
 b1 = 1._DP/6._DP; b2 = 1._DP/3._DP; b3 = 1._DP/3._DP; b4 = 1._DP/6._DP
 !--------------------------------------------------------------------!

 Q0 = Q
 ctime_n = ctime

 !Stage 1
 CALL compresid
 RK_k1 = resid
 Q = Q0 + 0.5_DP*dt*RK_k1
 ctime = ctime_n + 0.5_DP*dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 !Stage 2
 CALL compresid
 RK_k2 = resid
 Q = Q0 + 0.5_DP*dt*RK_k2

 !Stage 3
 CALL compresid
 RK_k3 = resid
 Q = Q0 + dt*RK_k3
 ctime = ctime_n + dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 !Stage 4
 CALL compresid
 RK_k4 = resid
 Q = Q0 + dt*(b1*RK_k1 + b2*RK_k2 + b3*RK_k3 + b4*RK_k4)
 
 END SUBROUTINE STEP_RK4
!----------------------------------------------------------------------------!


 SUBROUTINE STEP_SSPRK5

 USE setup
 USE MESH2D
 IMPLICIT NONE

 REAL(KIND=DP), DIMENSION(4,N,N,NCELL) :: Q0, Q1, Q2, Q3, Q4
 REAL(KIND=DP) :: a20, a21
 REAL(KIND=DP) :: a30, a31, a32
 REAL(KIND=DP) :: a40, a41, a42, a43
 REAL(KIND=DP) :: a50, a51, a52, a53, a54
 REAL(KIND=DP) :: c1, c2, c3, c4, c5, c54
 REAL(KIND=DP) :: ctime_n

 !-------------------------------------------------------!
 !Define RK coefficients:
 c1 = 0.391752226571890_DP
 c2 = 0.368410593050371_DP
 c3 = 0.251891774271694_DP
 c4 = 0.544974750228521_DP
 c5 = 0.226007483236906_DP; c54 = 0.063692468666290_DP

 a20 = 0.444370493651235_DP; a21 = 0.555629506348765_DP

 a30 = 0.620101851488403_DP; a31 = 0._DP
 a32 = 0.379898148511597_DP

 a40 = 0.178079954393132_DP; a41 = 0._DP; a42 = 0._DP
 a43 = 0.821920045606868_DP

 a50 = 0._DP; a51 = 0._DP; a52 = 0.517231671970585_DP
 a53 = 0.096059710526147_DP; a54 = 0.386708617503269_DP
 !-------------------------------------------------------!

 Q0 = Q;
 ctime_n = ctime

 !Stage 1
 CALL compresid
 RK_k1 = resid
 Q = Q0 + c1*dt*RK_k1
 Q1 = Q
 ctime = ctime_n + 0.391752226571890_DP*dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 !Stage 2
 CALL compresid
 RK_k2 = resid
 Q = a20*Q0 + a21*Q1 + c2*dt*RK_k2
 Q2 = Q
 ctime = ctime_n + 0.586079688967797_DP*dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 !Stage 3
 CALL compresid
 RK_k3 = resid
 Q = a30*Q0 + a31*Q1 + a32*Q2 + c3*dt*RK_k3
 Q3 = Q
 ctime = ctime_n + 0.474542363026867_DP*dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 !Stage 4
 CALL compresid
 RK_k4 = resid
 Q = a40*Q0 + a41*Q1 + a42*Q2 + a43*Q3 + c4*dt*RK_k4
 Q4 = Q
 ctime = ctime_n + 0.689388353336143_DP*dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS

 !Stage 5
 CALL compresid
 RK_k5 = resid
 Q = a50*Q0 + a51*Q1 + a52*Q2 + a53*Q3 + a54*Q4 + c54*RK_k4 + c5*RK_k5
 ctime = ctime_n + dt
 IF (IS_DEFORMING) CALL UPDATE_MESH_METRICS
 
 END SUBROUTINE STEP_SSPRK5
 !----------------------------------------------------------------------------!


 SUBROUTINE check_min_density

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER :: ic

 CHARACTER (LEN = 100) :: ERR_MSG

 DO ic = 1, NCELL
	IF ( MINVAL(Q(1,:,:,ic)) .LT. 0._DP ) THEN
		WRITE(ERR_MSG,'(A,I4)') ' NEGATIVE DENSITY DETECTED AT CELL: ', ic
		CALL SHUTDOWN (ERR_MSG)
	END IF
 END DO
 
 END SUBROUTINE check_min_density
 !----------------------------------------------------------------------------!


 SUBROUTINE comp_error

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER :: IC, is, js

 REAL(KIND=DP) :: x, y
 REAL(KIND=DP) :: r_in, x0, y0, rho_in, Mach_in, p_in
 REAL(KIND=DP) :: Xmin, Xmax, Ymin, Ymax
 REAL(KIND=DP) :: rhor, temp, r_vor, velinf, Machp
 REAL(KIND=DP) :: L1_error, L2_error, error


 L1_error = 0.0_DP; L2_error = 0.0_DP

 IF (INIT_COND .EQ. SS_VORTEX) THEN

	!Define problem parameters:
	r_in = 1.0; x0 = 0.0; y0 = 0.0
	rho_in = 1.0; Mach_in = 2.25

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				x = XXSolu(1,is,js,IC); y = XXSolu(2,is,js,IC)
				temp  = (r_in**2)/( (x-x0)**2+(y-y0)**2 )
				rhor  = rho_in*( 1.+0.5*(gam-1.0)*Mach_in**2*(1.-temp) )**(1.0/(gam-1.0))

				error = ABS(rhor - Q(1,is,js,IC))
				L1_error = L1_error + error
				L2_error = L2_error + error**2
			END DO
		END DO
	END DO

 ELSEIF (INIT_COND .EQ. X_SINWAVE) THEN

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				x = XXSolu(1,is,js,IC); y = XXSolu(2,is,js,IC)
				rhor = rhoinf + 0.1_DP*SIN( pi*((x-uinf*ctime)+(y-vinf*ctime)) )

				error = ABS(rhor - Q(1,is,js,IC))
				L1_error = L1_error + error
				L2_error = L2_error + error**2
			END DO
		END DO
	END DO

 ELSEIF (INIT_COND .EQ. EU_VORTEX) THEN

	!Define problem parameters
	x0 = uinf*ctime; y0 = vinf*ctime

	Machp = 0.3_DP						!Mach number of problem
	r_vor = 1._DP						!Radius of vortex
	velinf = SQRT(uinf**2 + vinf**2)	!Total free-stream velocity

	Xmin = -5._DP; Xmax = 5._DP
	Ymin = -5._DP; Ymax = 5._DP
 
	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				x = XXSolu(1,is,js,IC); y = XXSolu(2,is,js,IC)
				IF (x .LT. Xmin+x0) x = Xmax + x

				temp = ((x-x0)**2 + (y-y0)**2)/r_vor**2

				rhor   = rhoinf*( 1._DP-((gam-1._DP)*(Machp**2)/2._DP)* &
					   EXP(-1._DP*temp) )**(1.0/(gam-1._DP))

				error = ABS(rhor - Q(1,is,js,IC))
				L1_error = L1_error + error
				L2_error = L2_error + error**2
			END DO
		END DO
	END DO

 ELSEIF (INIT_COND .EQ. FREE_STRM) THEN
	
	rhor = 1._DP

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				error = ABS(rhor - Q(1,is,js,IC))
				L1_error = L1_error + error
				L2_error = L2_error + error**2
			END DO
		END DO
	END DO

 END IF

 L1_error = L1_error/(N*N*NCELL)
 L2_error = SQRT(L2_error/(N*N*NCELL))

 WRITE(*,'(A,E22.15)') ' L1 Err Norm: ', L1_error
 WRITE(*,'(A,E22.15)') ' L2 Err Norm: ', L2_error
 PRINT*,''

 END SUBROUTINE comp_error
 !----------------------------------------------------------------------------!


 SUBROUTINE resid_norm

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER :: IC, is, js

 REAL(KIND=DP) :: residnorm

 residnorm = 0._DP

 DO IC = 1, NCELL
	DO is = 1, N
		DO js = 1, N
			residnorm = residnorm + resid(1,is,js,IC)**2
		END DO
	END DO
 END DO

 residnorm = SQRT(residnorm)/(N*N*NCELL)

 WRITE(*,'(A,I5,A,F10.5)') ' >>ITERATION: ', ITER, ' TIME ', ctime 
 WRITE(*,'(A,E22.15)') ' RESID  Norm: ', residnorm


 END SUBROUTINE resid_norm
 !----------------------------------------------------------------------------!
