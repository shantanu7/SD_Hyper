 SUBROUTINE interface_flux

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: ifn, ifpl, jfpl, ifpr, jfpr, nfp, sp
 INTEGER :: ic, iface, icleft, icright, ifacelc, ifacerc
 INTEGER :: faml, famr, sign_l, sign_r

 REAL(KIND=DP), DIMENSION(4) :: Qs, Qfl, Qfr, Fcl, Fcr
 REAL(KIND=DP), DIMENSION(3) :: normf

 DO ifn = 1, NINTER

	iface = IBFINTER(ifn)

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
		IF (faml .EQ. 1) THEN !xi (vertical) face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpl,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(ifpl,sp)
			END DO
			normf(1:3) = S1(1,1:3,ifpl,jfpl,ic)
			IF (ifpl .EQ. 1) sign_l = -1 	!local face# = 4
		ELSE				 !eta (horizontal) face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpl,sp,ic)
				Qfl(1:4) = Qfl(1:4) + Qs(1:4)*Lmat(jfpl,sp)
			END DO
			normf(1:3) = S2(2,1:3,ifpl,jfpl,ic)
			IF (jfpl .EQ. 1) sign_l = -1 	!local face# = 1
		END IF

		!Construct Q at FP for icright
		ifpr = iface2fp(N-nfp+1,ifacerc)
		jfpr = jface2fp(N-nfp+1,ifacerc)

		ic = icright
		sign_r = 1

		Qfr = 0._DP
		IF (famr .EQ. 1) THEN 	!xi face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,sp,jfpr,ic)
				Qfr(1:4) = Qfr(1:4) + Qs(1:4)*Lmat(ifpr,sp)
			END DO
!			normf(1:3) = S1(1,1:3,ifpr,jfpr,ic)
			IF (ifpr .EQ. N+1) sign_r = -1	!local face# = 2
		ELSE					!eta face
			DO sp = 1, N
				Qs(1:4)  = Q(1:4,ifpr,sp,ic)
				Qfr(1:4) = Qfr(1:4) + Qs(1:4)*Lmat(jfpr,sp)
			END DO
!			normf(1:3) = S2(2,1:3,ifpr,jfpr,ic)
			IF (jfpr .EQ. N+1) sign_r = -1	!local face# = 3
		END IF

		CALL Rusanov_Flux(Qfl,Qfr,Fcl,Fcr,sign_l,sign_r,normf)

		IF (faml .EQ. 1) THEN
			F1(1:4,ifpl,jfpl,icleft) = Fcl
		ELSE
			G2(1:4,ifpl,jfpl,icleft) = Fcl
		END IF

		IF (famr .EQ. 1) THEN
			F1(1:4,ifpr,jfpr,icright) = Fcr
		ELSE
			G2(1:4,ifpr,jfpr,icright) = Fcr
		END IF

	END DO !Loop over nfp

 END DO !Loop over ifn

 END SUBROUTINE interface_flux
 !-------------------------------------------------------------------------------!


 SUBROUTINE Rusanov_Flux(Qfl,Qfr,Fcl,Fcr,sign_l,sign_r,normf)

 USE MESH2D
 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), DIMENSION(4), INTENT(IN) :: Qfl, Qfr
 REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: normf
 REAL(KIND=DP), DIMENSION(4) :: Fnl, Fnr, Fcl, Fcr
 REAL(KIND=DP) :: rho_l, u_l, v_l, p_l, Vn_l, aVn_l
 REAL(KIND=DP) :: rho_r, u_r, v_r, p_r, Vn_r, aVn_r
 REAL(KIND=DP) :: eigv, c_a, magnorm

 INTEGER :: sign_l, sign_r

 !Get primitive vars from left cell
 rho_l 	= Qfl(1)
 u_l	= Qfl(2)/rho_l
 v_l	= Qfl(3)/rho_l
 p_l	= (gam-1.0_DP)*(Qfl(4) - 0.5_DP*rho_l*(u_l**2 + v_l**2))
 Vn_l	= u_l*normf(1) + v_l*normf(2) + normf(3)
 aVn_l	= u_l*normf(1) + v_l*normf(2)

 !Get primitive vars from right cell
 rho_r	= Qfr(1)
 u_r	= Qfr(2)/rho_r
 v_r	= Qfr(3)/rho_r
 p_r	= (gam-1.0_DP)*(Qfr(4) - 0.5_DP*rho_r*(u_r**2 + v_r**2))
 Vn_r	= u_r*normf(1) + v_r*normf(2) + normf(3)
 aVn_r	= u_r*normf(1) + v_r*normf(2)

 !Compute left flux vector
 Fnl(1)	= Qfl(1)*Vn_l
 Fnl(2) = Qfl(2)*Vn_l + normf(1)*p_l
 Fnl(3) = Qfl(3)*Vn_l + normf(2)*p_l
 Fnl(4) = Qfl(4)*Vn_l + p_l*aVn_l

 !Compute right flux vector
 Fnr(1)	= Qfr(1)*Vn_r
 Fnr(2) = Qfr(2)*Vn_r + normf(1)*p_r
 Fnr(3) = Qfr(3)*Vn_r + normf(2)*p_r
 Fnr(4) = Qfr(4)*Vn_r + p_r*aVn_r

 !Common quantities
 magnorm= SQRT(normf(1)**2 + normf(2)**2)
 c_a	= SQRT( gam*(p_l+p_r)/(rho_l+rho_r) )
 eigv	= 0.5_DP*ABS(Vn_l + Vn_r) + c_a*magnorm
 
 !Compute interfacial Rusanov fluxes
 Fcl 	= 0.5_DP*(Fnl+Fnr - eigv*(Qfr-Qfl)*REAL(sign_l,KIND=DP))
 Fcr 	= 0.5_DP*(Fnl+Fnr - eigv*(Qfr-Qfl)*REAL(sign_r,KIND=DP))

 END SUBROUTINE Rusanov_Flux
 !-------------------------------------------------------------------------------!

