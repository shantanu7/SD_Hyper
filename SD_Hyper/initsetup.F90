 SUBROUTINE INITSETUP

 USE MESH2D
 USE setup
 IMPLICIT NONE

 CHARACTER (LEN = 100) :: ERR_MSG

 ALLOCATE(Q(4,N,N,NCELL), resid(4,N,N,NCELL))
 ALLOCATE(F1(4,N+1,N,NCELL), G2(4,N,N+1,NCELL))

 !Define positions of solution and flux points
 CALL SET_SOL_FLUX_POINTS

 CALL COMPUTE_L_M_MATRICES

 CALL FACE2FP_CONNEX

 CALL DISP_SP_FP_REPORT

 CALL MAPDER

 CALL INIT_MESH_METRICS

 CALL SET_INIT_COND

 IF (OUT_FORMAT .EQ. T_ASCII) THEN
	 CALL WRITE_TEC_DATA_NODAL
 ELSEIF (OUT_FORMAT .EQ. T_BINARY) THEN
	 CALL WRITE_TEC_DATA_NODAL_BINARY
 ELSE
	WRITE(ERR_MSG,'(A,I2)') ' INVALID OUTPUT FORMAT: ', OUT_FORMAT
	CALL SHUTDOWN (ERR_MSG)
 END IF

 END SUBROUTINE INITSETUP
 !-------------------------------------------------------------------------------!


 SUBROUTINE SET_SOL_FLUX_POINTS

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: isp, ifp

 !Allocate memory for solution and flux point structures
 ALLOCATE(Xs(N), Xf(N+1))

 !Define solution points -- Chebyshev-Gauss-Lobatto
 WRITE(*, '(A)', ADVANCE = 'NO') ' DEFINING SOLUTION POINTS... '

 DO isp = 1, N
	Xs(isp) = 0.5_DP*( 1._DP-COS( REAL(2._DP*isp-1._DP,KIND=16)/ &
								(2._DP*REAL(N,KIND=16))*pi ) )
 END DO

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''

 !Define position of flux points -- Legendre zeros + ends
 WRITE(*, '(A)', ADVANCE = 'NO') ' DEFINING FLUX POINTS... '

 Xf(1) = 0._DP; Xf(N+1) = 1._DP
 IF     (N .EQ. 2) THEN
	Xf(2) = 0.5_DP
 ELSEIF (N .EQ. 3) THEN
	Xf(2) = 0.5_DP*(1._DP - SQRT(3._DP)/3._DP)
	Xf(3) = 0.5_DP*(1._DP + SQRT(3._DP)/3._DP)
 ELSEIF (N .EQ. 4) THEN
	Xf(2) = 0.5_DP*(1._DP - SQRT(15._DP)/5._DP)
	Xf(3) = 0.5_DP
	Xf(4) = 0.5_DP*(1._DP + SQRT(15._DP)/5._DP)
 ELSEIF (N .EQ. 5) THEN
	Xf(2) = 0.5_DP*(1._DP - SQRT(525._DP+70._DP*SQRT(30._DP))/35._DP)
	Xf(3) = 0.5_DP*(1._DP - SQRT(525._DP-70._DP*SQRT(30._DP))/35._DP)
	Xf(4) = 0.5_DP*(1._DP + SQRT(525._DP-70._DP*SQRT(30._DP))/35._DP)
	Xf(5) = 0.5_DP*(1._DP + SQRT(525._DP+70._DP*SQRT(30._DP))/35._DP)
 ELSEIF (N .EQ. 6) THEN
	Xf(2) = (1. - SQRT(245._DP+14._DP*SQRT(70._DP))/21._DP)/2._DP
	Xf(3) = (1. - SQRT(245._DP-14._DP*SQRT(70._DP))/21._DP)/2._DP
	Xf(4) =  0.5_DP
	Xf(5) = (1. + SQRT(245._DP-14._DP*SQRT(70._DP))/21._DP)/2._DP
	Xf(6) = (1. + SQRT(245._DP+14._DP*SQRT(70._DP))/21._DP)/2._DP
 ELSEIF (N .EQ. 7) THEN
	Xf(2) = (1._DP - 0.932469514203_DP)/2._DP
	Xf(3) = (1._DP - 0.661209386466_DP)/2._DP
	Xf(4) = (1._DP - 0.238619186083_DP)/2._DP
	Xf(5) = (1._DP + 0.238619186083_DP)/2._DP
	Xf(6) = (1._DP + 0.661209386466_DP)/2._DP
	Xf(7) = (1._DP + 0.932469514203_DP)/2._DP
 END IF

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''

 END SUBROUTINE SET_SOL_FLUX_POINTS
 !-------------------------------------------------------------------------------!


 SUBROUTINE COMPUTE_L_M_MATRICES

 USE setup
 IMPLICIT NONE

 INTEGER :: isp, jsp, ifp, s, k
 
 REAL(KIND=DP) :: xval, summ, nume, deno

 ALLOCATE(Lmat(N+1,N), Mmat(N+1,N))

 WRITE(*, '(A)', ADVANCE = 'NO') ' CONSTRUCTING INTERPOLATION STRUCTURES... '

 !Lmat
 DO ifp = 1, N+1
	xval = Xf(ifp)
	DO isp = 1, N
		Lmat(ifp,isp) = hhval(isp, xval)
	END DO
 END DO

 !Mmat
 DO ifp = 1, N+1
	DO isp = 1, N
		xval = Xs(isp); summ = 0._DP

		DO k = 1, N+1
			IF (k .NE. ifp) THEN
				nume = 1._DP; deno = 1._DP
				DO s = 1, N+1
					IF (s .NE. ifp .AND. s .NE. k) nume = nume*(xval-Xf(s))
					
					IF (ifp .NE. s) deno = deno*(Xf(ifp)-Xf(s))
				END DO
				summ = summ + nume/deno
			END IF
		END DO

		Mmat(ifp,isp) = summ
	END DO
 END DO

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''
!----------------------------------------------------------!
 
 CONTAINS

 REAL(KIND=DP) FUNCTION hhval(i, xval)

 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), INTENT (IN) :: xval
 REAL(KIND=DP) :: nume, deno

 INTEGER, INTENT (IN) :: i
 INTEGER :: s

 nume = 1.0; deno = 1.0

 DO s = 1, N
	IF (s .NE. i) THEN
		nume = nume*(xval-Xs(s))
		deno = deno*(Xs(i)-Xs(s))
	END IF
 END DO 

 hhval = nume/deno
  
 END FUNCTION hhval
 
 END SUBROUTINE COMPUTE_L_M_MATRICES
 !-------------------------------------------------------------------------------!


 SUBROUTINE FACE2FP_CONNEX

 USE setup
 IMPLICIT NONE

 INTEGER :: nfp, iface

 WRITE(*,'(A)', ADVANCE = 'NO') ' CREATING FACE2FP CONNECTIVITY... '
 !Create connectivity for local face numbers and flux points on the faces.

 ALLOCATE(iface2fp(N,4), jface2fp(N,4))

 iface = 1
 DO nfp = 1, N
	iface2fp(nfp,iface) = nfp
	jface2fp(nfp,iface) = 1
 END DO
 
 iface = 2
 DO nfp = 1, N
	iface2fp(nfp,iface) = N+1
	jface2fp(nfp,iface) = nfp
 END DO

 iface = 3
 DO nfp = 1, N
	iface2fp(nfp,iface) = N-nfp+1
	jface2fp(nfp,iface) = N+1
 END DO

 iface = 4
 DO nfp = 1, N
	iface2fp(nfp,iface) = 1
	jface2fp(nfp,iface) = N-nfp+1
 END DO

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''

 END SUBROUTINE FACE2FP_CONNEX
 !-------------------------------------------------------------------------------!


 SUBROUTINE SET_INIT_COND

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: IC, is, js

 REAL(KIND=DP) :: xfp, yfp, xxf(2,4), rrf(2), x, y
 REAL(KIND=DP) :: r_in, x0, y0, rho_in, Mach_in, p_in
 REAL(KIND=DP) :: rhor, ur, vr, temp, plocal, Machp, V2
 REAL(KIND=DP) :: r_vor, velinf

 WRITE(*, '(A)', ADVANCE = 'NO') ' DEFINING INITIAL CONDITIONS... '

 ITER = 0; ctime = 0._DP

 IF (INIT_COND .EQ. SS_VORTEX) THEN

	!Define problem parameters:
	r_in = 1._DP; x0 = 0._DP; y0 = 0._DP
	rho_in = 1._DP; Mach_in = 2.25_DP; p_in = 1._DP/gam

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				!Define coordinates
				x = XXSolu(1,is,js,IC); y = XXSolu(2,is,js,IC)

				!Define local parameters
				temp  = (r_in**2)/( (x-x0)**2+(y-y0)**2 )
				rhor  = rho_in*( 1._DP+0.5_DP*(gam-1._DP)*Mach_in**2*&
						(1._DP-temp) )**(1._DP/(gam-1._DP))
				Machp = Mach_in*SQRT(temp)
				plocal = ((rhor/rho_in)**(gam))*p_in

				!Compute velocities
				vr = Machp*(x-x0)/SQRT((x-x0)**2+(y-y0)**2)
				ur = -1._DP*Machp*(y-y0)/SQRT((x-x0)**2+(y-y0)**2)

				!Assign solution variable values
				Q(1,is,js,IC) = rhor
				Q(2,is,js,IC) = rhor*ur
				Q(3,is,js,IC) = rhor*vr
				Q(4,is,js,IC) = plocal/(gam-1._DP) + 0.5_DP*rhor*(ur**2+vr**2)
			END DO
		END DO
	END DO

 ELSEIF (INIT_COND .EQ. X_SINWAVE) THEN 

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				x = XXSolu(1,is,js,IC); y = XXSolu(2,is,js,IC)

				rhor = rhoinf + 0.1_DP*SIN(pi*(x+y))
				ur = uinf; vr = vinf
				plocal = pinf

				Q(1,is,js,IC) = rhor
				Q(2,is,js,IC) = rhor*ur
				Q(3,is,js,IC) = rhor*vr
				Q(4,is,js,IC) = plocal/(gam-1._DP) + 0.5_DP*rhor*(ur**2+vr**2)
			END DO
		END DO
	END DO

 ELSEIF (INIT_COND .EQ. EU_VORTEX) THEN 

	!Define problem parameters
	x0 = 0._DP; y0 = 0._DP				!Eye of the vortex
	Machp = 0.3_DP						!Mach number of problem
	r_vor = 1._DP						!Radius of vortex
	velinf = SQRT(uinf**2 + vinf**2)	!Total free-stream velocity

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				x = XXSolu(1,is,js,IC); y = XXSolu(2,is,js,IC)
				temp = ((x-x0)**2 + (y-y0)**2)/r_vor**2

				rhor   = rhoinf*( 1._DP-((gam-1._DP)*(Machp**2)/2._DP)* &
					   EXP(-1._DP*temp) )**(1.0/(gam-1._DP))
				ur     = uinf - velinf*(y-y0)/r_vor*EXP(-1._DP*temp/2._DP)
				vr     = vinf + velinf*(x-x0)/r_vor*EXP(-1._DP*temp/2._DP)
				plocal = pinf*( 1._DP-((gam-1._DP)*(Machp**2)/2._DP)* &
					   EXP(-1._DP*temp) )**(gam/(gam-1._DP))

				Q(1,is,js,IC) = rhor
				Q(2,is,js,IC) = rhor*ur
				Q(3,is,js,IC) = rhor*vr
				Q(4,is,js,IC) = plocal/(gam-1._DP) + 0.5_DP*rhor*(ur**2+vr**2)
			END DO
		END DO
	END DO

 ELSEIF (INIT_COND .EQ. FREE_STRM) THEN

	DO IC = 1, NCELL
		DO is = 1, N
			DO js = 1, N
				rhor = 1.0; ur = 1.0; vr = 0.0
				plocal = 1.0
				Q(1,is,js,IC) = rhor
				Q(2,is,js,IC) = rhor*ur
				Q(3,is,js,IC) = rhor*vr
				Q(4,is,js,IC) = plocal/(gam-1._DP) + 0.5_DP*rhor*(ur**2+vr**2)
			END DO
		END DO
	END DO

 END IF

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''

 END SUBROUTINE SET_INIT_COND
 !-------------------------------------------------------------------------------!
