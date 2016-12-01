 SUBROUTINE MAPDER

 USE setup
 IMPLICIT NONE

 INTEGER :: is, js, ifp, jfp

 REAL(KIND=DP) :: xi, eta

 !Allocate associated memory
 ALLOCATE(dpdxs(4,3,N,N))
 ALLOCATE(dpdxf1(4,3,N+1,N))
 ALLOCATE(dpdxf2(4,3,N,N+1))

 WRITE(*, '(A)', ADVANCE = 'NO') ' COMPUTING DERIVATIVES OF SHAPE FUNCTIONS... '

 !Define derivatives at Solution points
 DO is = 1, N
	xi = Xs(is)
	DO js = 1, N
		eta = Xs(js)

		dpdxs(1,1,is,js) = eta-1._DP
		dpdxs(1,2,is,js) = xi-1._DP
		dpdxs(1,3,is,js) = (1._DP-xi)*(1._DP-eta)

		dpdxs(2,1,is,js) = 1._DP-eta
		dpdxs(2,2,is,js) = -xi
		dpdxs(2,3,is,js) = xi*(1._DP-eta)

		dpdxs(3,1,is,js) = eta
		dpdxs(3,2,is,js) = xi
		dpdxs(3,3,is,js) = (xi)*(eta)

		dpdxs(4,1,is,js) = -eta
		dpdxs(4,2,is,js) = 1._DP-xi
		dpdxs(4,3,is,js) = (1._DP-xi)*eta
	END DO
 END DO

 !Define derivatives at first family of flux points
 DO ifp = 1, N+1
	xi = Xf(ifp)
	DO js = 1, N
		eta = Xs(js)

		dpdxf1(1,1,ifp,js) = eta-1._DP
		dpdxf1(1,2,ifp,js) = xi-1._DP
		dpdxf1(1,3,ifp,js) = (1._DP-xi)*(1._DP-eta)

		dpdxf1(2,1,ifp,js) = 1._DP-eta
		dpdxf1(2,2,ifp,js) = -xi
		dpdxf1(2,3,ifp,js) = xi*(1._DP-eta)

		dpdxf1(3,1,ifp,js) = eta
		dpdxf1(3,2,ifp,js) = xi
		dpdxf1(3,3,ifp,js) = (xi)*(eta)

		dpdxf1(4,1,ifp,js) = -eta
		dpdxf1(4,2,ifp,js) = 1._DP-xi
		dpdxf1(4,3,ifp,js) = (1._DP-xi)*eta
	END DO
 END DO

 !Define derivatives at second family of flux points
 DO is = 1, N
	xi = Xs(is)
	DO jfp = 1, N+1
		eta = Xf(jfp)

		dpdxf2(1,1,is,jfp) = eta-1._DP
		dpdxf2(1,2,is,jfp) = xi-1._DP
		dpdxf2(1,3,is,jfp) = (1._DP-xi)*(1._DP-eta)

		dpdxf2(2,1,is,jfp) = 1._DP-eta
		dpdxf2(2,2,is,jfp) = -xi
		dpdxf2(2,3,is,jfp) = xi*(1._DP-eta)

		dpdxf2(3,1,is,jfp) = eta
		dpdxf2(3,2,is,jfp) = xi
		dpdxf2(3,3,is,jfp) = (xi)*(eta)

		dpdxf2(4,1,is,jfp) = -eta
		dpdxf2(4,2,is,jfp) = 1._DP-xi
		dpdxf2(4,3,is,jfp) = (1._DP-xi)*eta
	END DO
 END DO

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''
 

 END SUBROUTINE MAPDER
 !-------------------------------------------------------------------------------!
