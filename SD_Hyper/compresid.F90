 SUBROUTINE compresid

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: IC, is, js, rfp
 
 REAL(KIND=DP), DIMENSION(4) :: Qs, dF_xi, dG_eta
 REAL(KIND=DP) :: source1, source2

 CALL compflux

 !Compute flux derivatives at solution points
 DO IC = 1, NCELL

	!Transfer fluxes from FP to SP and add GCL-compliance terms
	DO is = 1, N
		DO js = 1, N
			Qs(1:4) = Q(1:4,is,js,IC)

			dF_xi   = 0._DP; dG_eta  = 0._DP
			source1 = 0._DP; source2 = 0._DP

			DO rfp = 1, N+1
				dF_xi(1:4)  = dF_xi(1:4)  + F1(1:4,rfp,js,IC)*Mmat(rfp,is)
				dG_eta(1:4) = dG_eta(1:4) + G2(1:4,is,rfp,IC)*Mmat(rfp,js)
				IF (IS_DEFORMING) THEN
					!Source terms for deforming meshes
					source1 = source1 + S1(1,3,rfp,js,IC)*Mmat(rfp,is)
					source2 = source2 + S2(2,3,is,rfp,IC)*Mmat(rfp,js)
				END IF
			END DO

			resid(1:4,is,js,IC) = -(dF_xi(1:4) + dG_eta(1:4)) + &
								  Qs(1:4)*(source1 + source2)

			resid(1:4,is,js,IC) = resid(1:4,is,js,IC)/Jac(is,js,IC)
		END DO
	END DO

 END DO
 
 END SUBROUTINE compresid
 !----------------------------------------------------------------------------!
 
