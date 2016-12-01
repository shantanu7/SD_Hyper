SUBROUTINE DISP_INIT_MESSAGE

 USE setup, ONLY: IS_RESTART

 IF (.NOT. IS_RESTART) THEN
	 PRINT*,''
	 PRINT*,'!---------------------------------------------------------------------!'
	 PRINT*,'!                          SD_HYPER (NEW RUN)                         !'
	 PRINT*,'!---------------------------------------------------------------------!'
	 PRINT*,''
 END IF

 END SUBROUTINE DISP_INIT_MESSAGE
!----------------------------------------------------------------------------!


 SUBROUTINE DISP_GEOM_REPORT

 USE MESH2D

 PRINT*,'--------------------------------------------------'
 PRINT*,'GEOMETRY AND MESH REPORT:'
 PRINT*,''
 PRINT*, 'PROJECT NAME: ', TRIM(PROJECT_NAME)
 WRITE(*,'(A,I5)')' NUMBER OF VERTICES: ', NVERT
 WRITE(*,'(A,I5)')' NUMBER OF CELLS: ', NCELL
 WRITE(*,'(A,I5)')' TOTAL NUMBER OF FACES DETECTED: ', NFACES
 PRINT*,''

 WRITE(*,'(A,I5)')' NUMBER OF CELL INTERFACES DETECTED: ', NINTER
 PRINT*,''

 WRITE(*,'(A)')' BOUNDARY DATA READ:'
 WRITE(*,'(I4,1X,A)') NBOUN, 'BOUNDARY FACES'
 WRITE(*,'(2X,A,I4,1X,A)') '->', NINLE, 'INLE FACES'
 WRITE(*,'(2X,A,I4,1X,A)') '->', NOUTL, 'OUTL FACES'
 WRITE(*,'(2X,A,I4,1X,A)') '->', NSYMM, 'SYMM FACES'
 WRITE(*,'(2X,A,I4,1X,A)') '->', NPERO, 'PERO FACES'
 WRITE(*,'(2X,A,I4,1X,A)') '->', NUKWN, 'UKWN FACES'
 PRINT*,'--------------------------------------------------'
 PRINT*,''

 END SUBROUTINE DISP_GEOM_REPORT
!----------------------------------------------------------------------------!


 SUBROUTINE DISP_SP_FP_REPORT

 USE setup
 IMPLICIT NONE

 INTEGER :: ifp

 PRINT*,'--------------------------------------------------'
 WRITE(*,'(A,I1)'),' SD SCHEME REPORT FOR SPATIAL ORDER: ', N
 PRINT*,''
 PRINT*,'SOLUTION POINT LOCATIONS: CHEBYSHEV-GAUSS-LOBATTO'
 WRITE(*,'(3(1X,E22.15))') Xs(1:N)
 PRINT*,''
 PRINT*,'FLUX POINT LOCATIONS    : ZEROS OF LEGENDRE-POLY.'
 WRITE(*,'(3(1X,E22.15))') Xf(1:N+1)
 PRINT*,''
 PRINT*,'INTERPOLATION MATRICES:'
 PRINT*,'Lmat:'
 DO ifp = 1, N+1
	WRITE(*,'(3(1X,E22.15))') Lmat(ifp,:)
 END DO 
 PRINT*,''
 PRINT*,'Mmat:'
 DO ifp = 1, N+1
	WRITE(*,'(3(1X,E22.15))') Mmat(ifp,:)
 END DO 
 PRINT*,'--------------------------------------------------'
 PRINT*,''

 END SUBROUTINE DISP_SP_FP_REPORT
 !-------------------------------------------------------------------------------!


 SUBROUTINE SHUTDOWN (ERR_MSG)
 
 CHARACTER (LEN = 100), INTENT (IN) :: ERR_MSG

 PRINT*,'PROCESS INTERRUPTED! '
 PRINT*,''
 WRITE(*,'(A,A)')' FATAL ERROR: ', TRIM(ERR_MSG)
 WRITE(*,'(A)')' ABORTING! '
 PRINT*,''

 STOP

 END SUBROUTINE SHUTDOWN
!----------------------------------------------------------------------------!
