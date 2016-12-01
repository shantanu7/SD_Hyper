 SUBROUTINE CALCULATE_JACOBIAN

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: IC, is, js, ifp, jfp, k

 REAL(KIND=DP) :: x_xi, x_eta, x_tau
 REAL(KIND=DP) :: y_xi, y_eta, y_tau

 !Define Jacobian and its determinant at solution points
 DO IC = 1, NCELL

	 DO is = 1, N
		DO js = 1, N
			x_xi = 0._DP; x_eta = 0._DP; x_tau = 0._DP
			y_xi = 0._DP; y_eta = 0._DP; y_tau = 0._DP
			DO k = 1, 4
				x_xi  = x_xi  + dpdxs(k,1,is,js)*XVMP(IVCELL(IC,k))
				x_eta = x_eta + dpdxs(k,2,is,js)*XVMP(IVCELL(IC,k))
				x_tau = x_tau + dpdxs(k,3,is,js)*XVMP_VEL(IVCELL(IC,k))

				y_xi  = y_xi  + dpdxs(k,1,is,js)*YVMP(IVCELL(IC,k))
				y_eta = y_eta + dpdxs(k,2,is,js)*YVMP(IVCELL(IC,k))
				y_tau = y_tau + dpdxs(k,3,is,js)*YVMP_VEL(IVCELL(IC,k))
			END DO

			Jac(is,js,IC) = x_xi*y_eta - x_eta*y_xi
		END DO
	 END DO

 	!Define adjoint of Jacobian at flux points (family 1)
	DO ifp = 1, N+1
		DO js = 1, N
			x_xi = 0._DP; x_eta = 0._DP; x_tau = 0._DP
			y_xi = 0._DP; y_eta = 0._DP; y_tau = 0._DP
			DO k = 1, 4
				x_xi  = x_xi  + dpdxf1(k,1,ifp,js)*XVMP(IVCELL(IC,k))
				x_eta = x_eta + dpdxf1(k,2,ifp,js)*XVMP(IVCELL(IC,k))
				x_tau = x_tau + dpdxf1(k,3,ifp,js)*XVMP_VEL(IVCELL(IC,k))

				y_xi  = y_xi  + dpdxf1(k,1,ifp,js)*YVMP(IVCELL(IC,k))
				y_eta = y_eta + dpdxf1(k,2,ifp,js)*YVMP(IVCELL(IC,k))
				y_tau = y_tau + dpdxf1(k,3,ifp,js)*YVMP_VEL(IVCELL(IC,k))
			END DO

			S1(1,1,ifp,js,IC) =  y_eta
			S1(1,2,ifp,js,IC) = -x_eta
			S1(1,3,ifp,js,IC) =  x_eta*y_tau - x_tau*y_eta

			S1(2,1,ifp,js,IC) = -y_xi
			S1(2,2,ifp,js,IC) =  x_xi
			S1(2,3,ifp,js,IC) =  x_tau*y_xi - x_xi*y_tau
		END DO
	END DO

	!Define adjoint Jacobian at flux points (family 2)
	DO is = 1, N
		DO jfp = 1, N+1
			x_xi = 0._DP; x_eta = 0._DP; x_tau = 0._DP
			y_xi = 0._DP; y_eta = 0._DP; y_tau = 0._DP
			DO k = 1, 4
				x_xi  = x_xi  + dpdxf2(k,1,is,jfp)*XVMP(IVCELL(IC,k))
				x_eta = x_eta + dpdxf2(k,2,is,jfp)*XVMP(IVCELL(IC,k))
				x_tau = x_tau + dpdxf2(k,3,is,jfp)*XVMP_VEL(IVCELL(IC,k))

				y_xi  = y_xi  + dpdxf2(k,1,is,jfp)*YVMP(IVCELL(IC,k))
				y_eta = y_eta + dpdxf2(k,2,is,jfp)*YVMP(IVCELL(IC,k))
				y_tau = y_tau + dpdxf2(k,3,is,jfp)*YVMP_VEL(IVCELL(IC,k))
			END DO

			S2(1,1,is,jfp,IC) =  y_eta
			S2(1,2,is,jfp,IC) = -x_eta
			S2(1,3,is,jfp,IC) =  x_eta*y_tau - x_tau*y_eta

			S2(2,1,is,jfp,IC) = -y_xi
			S2(2,2,is,jfp,IC) =  x_xi
			S2(2,3,is,jfp,IC) =  x_tau*y_xi - x_xi*y_tau
		END DO
	END DO

 END DO

 END SUBROUTINE CALCULATE_JACOBIAN
