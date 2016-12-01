 SUBROUTINE BC_flux

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: ifn, ifpl, jfpl, ifpr, jfpr, nfp, sp
 INTEGER :: ic, iface, icleft, icright, ifacelc, ifacerc
 INTEGER :: faml, famr, sign_l, sign_r
 INTEGER :: k

 REAL(KIND=DP), DIMENSION(4) :: Qs, Qfl, Qfr, Fcl, Fcr
 REAL(KIND=DP), DIMENSION(3) :: normf
 REAL(KIND=DP) :: magnorm
 REAL(KIND=DP) :: rhow, uw, vw, VN
 REAL(KIND=DP) :: xfp, yfp, xxf(2,4), rrf(2), x, y
 REAL(KIND=DP) :: r_in, x0, y0, rho_in, Mach_in, p_in
 REAL(KIND=DP) :: rhor, ur, vr, temp, plocal, Machp, V2

 !Inlet faces:
 !-----------
 DO ifn = 1, NINLE

	iface   = IBFINLE(ifn)
	icleft  = IF2C(iface,1)
	ifacelc = IF2C(iface,3)

	faml = MOD(ifacelc,2) + 1

	!Get coordinates of corner vertices.
	DO k = 1, 4
		xxf(1,k) = XVMP(IVCELL(icleft,k))
		xxf(2,k) = YVMP(IVCELL(icleft,k))
	END DO

	DO nfp = 1, N
		
		!Construct Q at FP for icleft
		ifpl = iface2fp(nfp,ifacelc)
		jfpl = jface2fp(nfp,ifacelc)

		ic = icleft
		sign_l = 1
		
		Qfl = 0._DP
		IF (faml .EQ. 1) THEN	!xi-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpl,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(ifpl,sp)
			END DO
			normf(1:3) = S1(1,1:3,ifpl,jfpl,ic)
			IF (ifpl .EQ. 1) sign_l = -1 
		ELSE					!eta-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpl,sp,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(jfpl,sp)
			END DO
			normf(1:3) = S2(2,1:3,ifpl,jfpl,ic)
			IF (jfpl .EQ. 1) sign_l = -1
		END IF

		magnorm = SQRT(normf(1)**2 + normf(2)**2)

		IF (INIT_COND .EQ. SS_VORTEX) THEN

			!Define problem parameters:
			r_in = 1.0_DP; x0 = 0.0_DP; y0 = 0.0_DP
			rho_in = 1.0_DP; Mach_in = 2.25_DP; p_in = 1.0_DP/gam

			IF (faml .EQ. 1) THEN
				xfp = Xf(ifpl); yfp = Xs(jfpl)
			ELSE
				xfp = Xs(ifpl); yfp = Xf(jfpl)
			END IF

			CALL XYCoor_AT_GPS(4,xxf,rrf,xfp,yfp)

			x = rrf(1); y = rrf(2)

			!Define flow parameters
			temp  = (r_in**2)/( (x-x0)**2+(y-y0)**2 )
			rhor  = rho_in*( 1._DP+0.5_DP*(gam-1._DP)*Mach_in**2* & 
					(1._DP-temp) )**(1._DP/(gam-1._DP))
			Machp = Mach_in*SQRT(temp)
			plocal = ((rhor/rho_in)**gam)*p_in
		   
			!Compute velocities
			vr = 		Machp*(x-x0)/SQRT((x-x0)**2+(y-y0)**2)
			ur = -1._DP*Machp*(y-y0)/SQRT((x-x0)**2+(y-y0)**2)

			!Define right-side solution vector
			Qfr(1) = rhor
			Qfr(2) = rhor*ur
			Qfr(3) = rhor*vr
			Qfr(4) = plocal/(gam-1._DP) + 0.5_DP*rhor*(ur**2+vr**2)

		END IF

		sign_r = sign_l

		CALL Rusanov_Flux(Qfl,Qfr,Fcl,Fcr,sign_l,sign_r,normf)
		
		IF (faml .EQ. 1) THEN
			F1(1:4,ifpl,jfpl,icleft) = Fcl
		ELSE
			G2(1:4,ifpl,jfpl,icleft) = Fcl
		END IF

	END DO !Loop over nfp

 END DO !Loop over ifn


 !Outlet faces:
 !------------
 DO ifn = 1, NOUTL

	iface   = IBFOUTL(ifn)
	icleft  = IF2C(iface,1)
	ifacelc = IF2C(iface,3)

	faml = MOD(ifacelc,2) + 1

	DO nfp = 1, N
		
		!Construct Q at FP for icleft
		ifpl = iface2fp(nfp,ifacelc)
		jfpl = jface2fp(nfp,ifacelc)

		ic = icleft
		sign_l = 1; sign_r = 1

		Qfl = 0._DP
		IF (faml .EQ. 1) THEN	!xi-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpl,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(ifpl,sp)
			END DO
			normf(1:3) = S1(1,1:3,ifpl,jfpl,ic)
			IF (ifpl .EQ. 1) sign_l = -1 
		ELSE					!eta-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpl,sp,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(jfpl,sp)
			END DO
			normf(1:3) = S2(2,1:3,ifpl,jfpl,ic)
			IF (jfpl .EQ. 1) sign_l = -1
		END IF

		magnorm = SQRT(normf(1)**2 + normf(2)**2)

		sign_r = sign_l

		!Define right-side solution vector
		!Extrapolation outlet BC -- supersonic flows
		Qfr(1) = Qfl(1)
		Qfr(2) = Qfl(2)
		Qfr(3) = Qfl(3)
		Qfr(4) = Qfl(4)

		CALL Rusanov_Flux(Qfl,Qfr,Fcl,Fcr,sign_l,sign_r,normf)

		IF (faml .EQ. 1) THEN
			F1(1:4,ifpl,jfpl,icleft) = Fcl
		ELSE 
			G2(1:4,ifpl,jfpl,icleft) = Fcl
		END IF

	END DO

 END DO


 !Symmetric faces:
 !---------------
 DO ifn = 1, NSYMM
	
	iface   = IBFSYMM(ifn)
	icleft  = IF2C(iface,1)
	ifacelc = IF2C(iface,3)

	faml = MOD(ifacelc,2) + 1

	DO nfp = 1, N
		
		!Construct Q at FP for icleft
		ifpl = iface2fp(nfp,ifacelc)
		jfpl = jface2fp(nfp,ifacelc)

		ic = icleft
		sign_l = 1; sign_r = 1
		
		Qfl = 0._DP
		IF (faml .EQ. 1) THEN	!xi-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpl,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(ifpl,sp)
			END DO
			normf(1:3) = S1(1,1:3,ifpl,jfpl,ic)
			IF (ifpl .EQ. 1) sign_l = -1 
		ELSE					!eta-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpl,sp,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(jfpl,sp)
			END DO
			normf(1:3) = S2(2,1:3,ifpl,jfpl,ic)
			IF (jfpl .EQ. 1) sign_l = -1
		END IF

		magnorm = SQRT(normf(1)**2 + normf(2)**2)

		!Define Qfr
		rhow = Qfl(1)
		uw = Qfl(2)/rhow
		vw = Qfl(3)/rhow
		VN = (uw*normf(1) + vw*normf(2))/magnorm

		!Define right-side solution vector		
		Qfr(1) = Qfl(1)
		Qfr(2) = rhow*(uw - 2.0*VN*normf(1)/magnorm)
		Qfr(3) = rhow*(vw - 2.0*VN*normf(2)/magnorm)
		Qfr(4) = Qfl(4)

		CALL Rusanov_Flux(Qfl,Qfr,Fcl,Fcr,sign_l,sign_r,normf)
		
		IF (faml .EQ. 1) THEN
			F1(1:4,ifpl,jfpl,icleft) = Fcl
		ELSE 
			G2(1:4,ifpl,jfpl,icleft) = Fcl
		END IF

	END DO

 END DO


 !Periodic faces:
 !--------------
 DO ifn = 1, NPERO

	iface   = IBFPERO(ifn)

	icleft  = IF2C(iface,1); icright = IF2C(iface,2)
	ifacelc = IF2C(iface,3); ifacerc = IF2C(iface,4)

	faml = MOD(ifacelc,2) + 1; famr = MOD(ifacerc,2) + 1

	DO nfp = 1, N
		
		!Construct Q at FP for icleft
		ifpl = iface2fp(nfp,ifacelc)
		jfpl = jface2fp(nfp,ifacelc)

		ic = icleft
		sign_l = 1

		Qfl = 0._DP
		IF (faml .EQ. 1) THEN	!xi-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpl,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(ifpl,sp)
			END DO
			normf(1:3) = S1(1,1:3,ifpl,jfpl,ic)
			IF (ifpl .EQ. 1) sign_l = -1 
		ELSE					!eta-face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpl,sp,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(jfpl,sp)
			END DO
			normf(1:3) = S2(2,1:3,ifpl,jfpl,ic)
			IF (jfpl .EQ. 1) sign_l = -1
		END IF

		!Construct Q at FP for icright
		ifpr = iface2fp(N-nfp+1,ifacerc)
		jfpr = jface2fp(N-nfp+1,ifacerc)

		ic = icright
		sign_r = 1

		Qfr = 0._DP
		IF (famr .EQ. 1) THEN
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpr,ic)
				Qfr(1:4) = Qfr(1:4) + Qs(1:4)*Lmat(ifpr,sp)
			END DO
!			normf(1:3) = S1(1,1:3,ifpr,jfpr,ic)
			IF (ifpr .EQ. N+1) sign_r = -1
		ELSE
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpr,sp,ic)
				Qfr(1:4) = Qfr(1:4) + Qs(1:4)*Lmat(jfpr,sp)
			END DO
!			normf(1:3) = S2(2,1:3,ifpr,jfpr,ic)
			IF (jfpr .EQ. N+1) sign_r = -1 
		END IF

		CALL Rusanov_Flux(Qfl,Qfr,Fcl,Fcr,sign_l,sign_r,normf)

		IF (faml .EQ. 1) THEN
			F1(1:4,ifpl,jfpl,icleft) = Fcl
		ELSE 
			G2(1:4,ifpl,jfpl,icleft) = Fcl
		END IF

	END DO !Loop over nfp
	
 END DO !Loop over ifn

 
 END SUBROUTINE BC_flux
 !----------------------------------------------------------------------------!
