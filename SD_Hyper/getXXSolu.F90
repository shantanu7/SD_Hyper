 SUBROUTINE GetXXSolu

 USE MESH2D
 USE setup
 IMPLICIT NONE

 REAL(KIND=DP) :: xxf(2,4), xx(2)

 INTEGER :: IC, IV, k, is, js

 DO IC = 1, NCELL
	DO k = 1, 4
		xxf(1,k) = XVMP(IVCELL(IC,k))
		xxf(2,k) = YVMP(IVCELL(IC,k))
	END DO

	DO is = 1, N
		DO js = 1, N
			CALL XYCoor_AT_GPS(4,xxf,xx,Xs(is),Xs(js))
			XXSolu(1,is,js,IC) = xx(1)
			XXSolu(2,is,js,IC) = xx(2)
		END DO
	END DO
 END DO

 END SUBROUTINE GetXXSolu
 !----------------------------------------------------------------------------!


 SUBROUTINE XYCoor_AT_GPS(NV,xxf,xx,xi,eta)

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER, INTENT(IN) :: NV
 INTEGER :: i, iv
 
 REAL(KIND=DP) :: xx(2), xi, eta
 REAL(KIND=DP) :: Xis, Xjs, xxf(2,NV), func(NV)

 CALL MAP_BASE_FUNC(func,xi,eta)

 DO i = 1, 2
	xx(i) = 0.0
	DO iv = 1, NV
		xx(i) = xx(i) + xxf(i,iv)*func(iv)
	END DO
 END DO
 
 END SUBROUTINE XYCoor_AT_GPS
!----------------------------------------------------------------------------!


 SUBROUTINE MAP_BASE_FUNC(func,xi,eta)

 USE MESH2D
 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), INTENT(IN) :: xi, eta
 REAL(KIND=DP) :: func(4)

 func(1) = (1.0-xi)*(1.0-eta)
 func(2) = xi*(1.0-eta)
 func(3) = xi*eta
 func(4) = (1.0-xi)*eta
 
 END SUBROUTINE MAP_BASE_FUNC
!----------------------------------------------------------------------------!
