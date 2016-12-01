 SUBROUTINE compflux

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: IC, ifp, jfp, sp

 REAL(KIND=DP), DIMENSION(4):: Qs, Qf, F, G

 DO IC = 1, NCELL

	!Reconstruct fluxes at interior flux points (Family 1)
	DO ifp = 2, N
		DO jfp = 1, N
			Qf = 0._DP
			!Interpolate conserved variables from sol points.
			DO sp = 1, N
				Qs(1:4) = Q(1:4,sp,jfp,IC)
				Qf(1:4) = Qf(1:4) + Qs(1:4)*Lmat(ifp,sp)
			END DO
			!Compute inviscid fluxes in physical domain
			CALL getIfluxVectors(Qf,F,G)
			!Compute inviscid fluxes in computational domain
			F1(1:4,ifp,jfp,IC) = S1(1,1,ifp,jfp,IC)*F(1:4) + &
								 S1(1,2,ifp,jfp,IC)*G(1:4) + &
								 S1(1,3,ifp,jfp,IC)*Qf(1:4)
		END DO
	END DO

 	!Reconstruct fluxes at interior flux points (Family 2)
	DO ifp = 1, N
		DO jfp = 2, N
			Qf = 0._DP
			!Interpolate conserved variables from sol points
			DO sp = 1, N
				Qs(1:4) = Q(1:4,ifp,sp,IC)
				Qf(1:4) = Qf(1:4) + Qs(1:4)*Lmat(jfp,sp)
			END DO
			!Compute inviscid fluxes in physical domain
			CALL getIfluxVectors(Qf,F,G)
			!Compute inviscid fluxes in computational domain
			G2(1:4,ifp,jfp,IC) = S2(2,1,ifp,jfp,IC)*F(1:4) + &
								 S2(2,2,ifp,jfp,IC)*G(1:4) + &
								 S2(2,3,ifp,jfp,IC)*Qf(1:4)
		END DO
	END DO

 END DO

 !Compute fluxes at non-boundary cell interfaces
 CALL interface_flux

 !Compute fluxes at boundary faces
 CALL BC_flux

 END SUBROUTINE compflux
 !----------------------------------------------------------------------------!


 SUBROUTINE getIfluxVectors (Qi, Fi, Gi)

 USE MESH2D
 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), INTENT (IN), DIMENSION(4) :: Qi
 REAL(KIND=DP), DIMENSION(4) :: Fi, Gi
 REAL(KIND=DP) :: rho, u, v, ene, p
 
 !Calculate primitive variables
 rho = Qi(1)
 u 	 = Qi(2)/Qi(1)
 v 	 = Qi(3)/Qi(1)
 ene = Qi(4)/Qi(1)

 p = (gam-1.0)*( Qi(4) - 0.5*(Qi(2)**2+Qi(3)**2)/Qi(1) )

 !Calculate flux vectors
 Fi(1) = rho*u
 Fi(2) = rho*u**2 + p
 Fi(3) = rho*u*v
 Fi(4) = (rho*ene + p)*u

 Gi(1) = rho*v
 Gi(2) = rho*u*v
 Gi(3) = rho*v**2 + p
 Gi(4) = (rho*ene + p)*v
 
 END SUBROUTINE getIfluxVectors
 !----------------------------------------------------------------------------!
