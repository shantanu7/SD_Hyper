 SUBROUTINE CONNECTIVITY

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER :: IC, ICC, IBF, IV, IFV1, IFV2, JFV1, JFV2
 INTEGER :: ICLV1, ICLV2, ICRV1, ICRV2
 INTEGER :: IV1, IV2, IV3, IV4
 INTEGER :: iface, jface, kface, BFCOUNT, DEFICIT
 INTEGER :: icleft, icright
 INTEGER :: J, K
 INTEGER :: IF2V_TEMP(MAXFACES,2), ICVFILL(NVERT), ICVERT(NVERT,6)
 INTEGER :: IBFIN_TEMP(MAXBOUN), IBFOU_TEMP(MAXBOUN), &
			IBFSY_TEMP(MAXBOUN), IBFPE_TEMP(MAXBOUN), &
			IBFNT_TEMP(MAXFACES)

 LOGICAL :: IS_UNIQUE, IS_BOUNDARY

 CHARACTER (LEN = 100) :: ERR_MSG

 WRITE(*, '(A)', ADVANCE = 'NO') ' CREATING FACE DATA... '

 !Define IVFACE
 IVFACE(1,1) = 1; IVFACE(1,2) = 2
 IVFACE(2,1) = 2; IVFACE(2,2) = 3
 IVFACE(3,1) = 3; IVFACE(3,2) = 4
 IVFACE(4,1) = 4; IVFACE(4,2) = 1

 NFACES = 0

 ALLOCATE (IC2F(NCELL,4))

 DO IC = 1, NCELL

	!Create IF2V structures
	DO iface = 1, 4
		IF (NFACES .LT. 4) THEN !First four faces are new, by default
			IS_UNIQUE = .TRUE.
		ELSE !Searching for new faces required
			IS_UNIQUE = .TRUE.
			DO jface = 1, NFACES
				IFV1 = IVCELL(IC,IVFACE(iface,1))
				IFV2 = IVCELL(IC,IVFACE(iface,2))
				JFV1 = IF2V_TEMP(jface,1)
				JFV2 = IF2V_TEMP(jface,2)
				IF ((IFV1 .EQ. JFV1 .AND. IFV2 .EQ. JFV2) .OR. &
					(IFV1 .EQ. JFV2 .AND. IFV2 .EQ. JFV1)) THEN !New face not detected
					IS_UNIQUE = .FALSE.
					IC2F(IC,iface) = jface
					EXIT
				END IF
			END DO
		END IF
		IF (IS_UNIQUE) THEN
			NFACES = NFACES + 1
			IF2V_TEMP(NFACES,1) = IVCELL(IC,IVFACE(iface,1))
			IF2V_TEMP(NFACES,2) = IVCELL(IC,IVFACE(iface,2))
			IC2F(IC,iface) = NFACES
		END IF
	END DO !iface
	
 END DO !IC

 ALLOCATE(IF2V(NFACES,2))
 ALLOCATE(CHARACTER(LEN=16) :: FACEMARKER(NFACES))

 IF2V(:,:) = IF2V_TEMP(:NFACES,:)

 !Test for Boundary Faces
 BFCOUNT = 0; 
 NINLE = 0; NOUTL = 0; NSYMM = 0; NPERO = 0
 NINTER = 0; 
 NUKWN = 0
 FACEMARKER = ''

 DO iface = 1, NFACES
	IS_BOUNDARY = .FALSE.
	IFV1 = IF2V(iface,1); IFV2 = IF2V(iface,2)
	DO IBF = 1, NBOUN
		JFV1 = IVBOUN(IBF,1); JFV2 = IVBOUN(IBF,2)

		IF (iface .EQ. IBF2F(IBF)) THEN  !Boundary face detected
			IS_BOUNDARY = .TRUE.
			BFCOUNT = BFCOUNT + 1
			FACEMARKER(iface) = TRIM(FACEMARKER_BOUN(IBF))

			!Update Boundary type counters:
			IF 	   (TRIM(FACEMARKER(iface)) .EQ. 'INLE') THEN
				NINLE = NINLE + 1
				IBFIN_TEMP(NINLE) = iface
			ELSEIF (TRIM(FACEMARKER(iface)) .EQ. 'OUTL') THEN
				NOUTL = NOUTL + 1
				IBFOU_TEMP(NOUTL) = iface
			ELSEIF (TRIM(FACEMARKER(iface)) .EQ. 'SYMM') THEN
				NSYMM = NSYMM + 1
				IBFSY_TEMP(NSYMM) = iface
			ELSEIF (TRIM(FACEMARKER(iface)) .EQ. 'PERO') THEN
				NPERO = NPERO + 1
				IBFPE_TEMP(NPERO) = iface
			ELSE
				NUKWN = NUKWN + 1
			END IF

			EXIT
		END IF
	END DO !IBF	

	IF (.NOT. IS_BOUNDARY) THEN
		NINTER = NINTER + 1
		FACEMARKER(iface) = 'INTER'
		IBFNT_TEMP(NINTER) = iface
	END IF
 END DO !iface

 !Error check for bad boundary faces
 IF (NUKWN .GT. 0) THEN
	WRITE(ERR_MSG,'(I4,A)') NUKWN,' FACE(S) DECTECTED WITH UNRECOGNIZED BC TYPE!'
	CALL SHUTDOWN (ERR_MSG)
 END IF

 !Error check for incomplete boundary face detection
 DEFICIT = NBOUN - (NINLE+NOUTL+NSYMM+NPERO)
 IF (DEFICIT .NE. 0) THEN
	WRITE(ERR_MSG,'(A,I4,A)') ' BOUNDARY TYPES NOT DETECTED FOR ', DEFICIT,' FACE(S)!'
	CALL SHUTDOWN (ERR_MSG)
 END IF

 !Assign global face numbers to boundary face indices
 ALLOCATE(IBFINLE(NINLE), IBFOUTL(NOUTL), IBFSYMM(NSYMM), IBFPERO(NPERO))
 ALLOCATE(IBFINTER(NINTER))
 IBFINLE(:)  = IBFIN_TEMP(:NINLE)
 IBFOUTL(:)  = IBFOU_TEMP(:NOUTL)
 IBFSYMM(:)  = IBFSY_TEMP(:NSYMM)
 IBFPERO(:)  = IBFPE_TEMP(:NPERO)
 IBFINTER(:) = IBFNT_TEMP(:NINTER)

 !Find cells to the left and right of face.
 ALLOCATE(IF2C(NFACES,4))
 IF2C = 0; ICVFILL = 0

 DO IC = 1, NCELL
	DO K = 1, 4
		IV = IVCELL(IC,K)
		ICVFILL(IV) = ICVFILL(IV) + 1
		ICVERT(IV,ICVFILL(IV)) = IC
	END DO
 END DO

! DO IV = 1, NVERT
!	WRITE(765,*) IV, ICVERT(IV,:ICVFILL(IV))
! END DO

 !Find icleft and icright for each face
 DO IC = 1, NCELL
	DO kface = 1, 4
		iface = IC2F(IC,kface)

		IF ( IF2C(iface, 1) .EQ. 0 ) THEN !Face has not been tested yet

			ICLV1 = IVCELL(IC,IVFACE(kface,1))
			ICLV2 = IVCELL(IC,IVFACE(kface,2))
			IF2C(iface,1) = IC

			DO J = 1, ICVFILL(ICLV1) !Loop over cells surrounding ICLV1
				ICC = ICVERT(ICLV1,J)
				IF (IC .NE. ICC) THEN
					DO K = 1, 4
						ICRV2 = IVCELL(ICC,K)
						IF (ICRV2 .EQ. ICLV2) THEN
							IF2C(iface,2) = ICC
							EXIT
						END IF
					END DO
				END IF
			END DO

		END IF

	END DO
 END DO

 !Set ifacelc and ifacerc
 DO iface = 1, NFACES
	icleft = IF2C(iface,1)
	ICLV1 = IF2V(iface,1); ICLV2 = IF2V(iface,2)

	IV1 = IVCELL(icleft,1); IV2 = IVCELL(icleft,2);
	IV3 = IVCELL(icleft,3); IV4 = IVCELL(icleft,4);

	IF 		(((IV1 .EQ. ICLV1) .AND. (IV2 .EQ. ICLV2)) .OR. &
			 ((IV1 .EQ. ICLV2) .AND. (IV2 .EQ. ICLV1))) THEN
			IF2C(iface,3) = 1
	ELSEIF 	(((IV2 .EQ. ICLV1) .AND. (IV3 .EQ. ICLV2)) .OR. &
			 ((IV2 .EQ. ICLV2) .AND. (IV3 .EQ. ICLV1))) THEN
			IF2C(iface,3) = 2
	ELSEIF 	(((IV3 .EQ. ICLV1) .AND. (IV4 .EQ. ICLV2)) .OR. &
			 ((IV3 .EQ. ICLV2) .AND. (IV4 .EQ. ICLV1))) THEN
			IF2C(iface,3) = 3
	ELSEIF 	(((IV1 .EQ. ICLV1) .AND. (IV4 .EQ. ICLV2)) .OR. &
			 ((IV1 .EQ. ICLV2) .AND. (IV4 .EQ. ICLV1))) THEN
			IF2C(iface,3) = 4
	END IF

	IF (IF2C(iface,2) .NE. 0) THEN
		icright = IF2C(iface,2)
		ICRV1 = IF2V(iface,1); ICRV2 = IF2V(iface,2)

		IV1 = IVCELL(icright,1); IV2 = IVCELL(icright,2);
		IV3 = IVCELL(icright,3); IV4 = IVCELL(icright,4);

		IF 		(((IV1 .EQ. ICRV1) .AND. (IV2 .EQ. ICRV2)) .OR. &
				 ((IV1 .EQ. ICRV2) .AND. (IV2 .EQ. ICRV1))) THEN
				IF2C(iface,4) = 1
		ELSEIF 	(((IV2 .EQ. ICRV1) .AND. (IV3 .EQ. ICRV2)) .OR. &
				 ((IV2 .EQ. ICRV2) .AND. (IV3 .EQ. ICRV1))) THEN
				IF2C(iface,4) = 2
		ELSEIF 	(((IV3 .EQ. ICRV1) .AND. (IV4 .EQ. ICRV2)) .OR. &
				 ((IV3 .EQ. ICRV2) .AND. (IV4 .EQ. ICRV1))) THEN
				IF2C(iface,4) = 3
		ELSEIF 	(((IV1 .EQ. ICRV1) .AND. (IV4 .EQ. ICRV2)) .OR. &
				 ((IV1 .EQ. ICRV2) .AND. (IV4 .EQ. ICRV1))) THEN
				IF2C(iface,4) = 4
		END IF
	END IF
 END DO

 !Diagnostic to test IF2C
! DO iface = 1, NFACES
!	WRITE(765,'(7(I3,2X))') iface, IF2V(iface,:), IF2C(iface,:4)
! END DO

 IF (NPERO .GT. 0) CALL MATCH_PERO_FACES

!Diagnostic to test IF2C + match_pero
! DO iface = 1, NFACES
!	WRITE(765,'(5(I3,2X))') iface, IF2V(iface,:), IF2C(iface,:2)
! END DO

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''
 
 END SUBROUTINE CONNECTIVITY
!----------------------------------------------------------------------------!


 SUBROUTINE MATCH_PERO_FACES

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER :: IBF, JBF, iface, jface
 INTEGER :: icleft, icright, ifacelc, ifacerc
 INTEGER :: jcleft, jcright
 INTEGER :: IV1, IV2, JV1, JV2

 REAL(KIND=DP), PARAMETER :: tol = 1.0E-3
 REAL(KIND=DP) :: err11, err22, err12, err21
 
 LOGICAL :: FOUND_MATCH

 CHARACTER (LEN = 100) :: ERR_MSG

 DO IBF = 1, NPERO
	iface = IBFPERO(IBF)
	IV1 = IF2V(iface,1); IV2 = IF2V(iface,2)

	IF (IF2C(iface,2) .EQ. 0) THEN !Do only if face not previously searched
		FOUND_MATCH = .FALSE.
		DO JBF = 1, NPERO
			IF (IBF .NE. JBF) THEN
				jface = IBFPERO(JBF)
				JV1 = IF2V(jface,1); JV2 = IF2V(jface,2)
				IF (YV(IV1) .EQ. YV(IV2)) THEN	!Horizontal face
					err11 = ABS(XV(IV1) - XV(JV1))
					err22 = ABS(XV(IV2) - XV(JV2))
					err12 = ABS(XV(IV1) - XV(JV2))
					err21 = ABS(XV(IV2) - XV(JV1))
				ELSE	!Vertical face
					err11 = ABS(YV(IV1) - YV(JV1))
					err22 = ABS(YV(IV2) - YV(JV2))
					err12 = ABS(YV(IV1) - YV(JV2))
					err21 = ABS(YV(IV2) - YV(JV1))
				END IF
				IF ( (err11.LT.tol .AND. err22.LT.tol) .OR. &
					 (err21.LT.tol .AND. err12.LT.tol) ) THEN
					FOUND_MATCH = .TRUE.
					IF2C(iface,2) = IF2C(jface,1)
					IF2C(iface,4) = IF2C(jface,3)
					!Match the j-face
					IF2C(jface,2) = IF2C(iface,1)
					IF2C(jface,4) = IF2C(iface,3)
					EXIT
				END IF
			END IF
		END DO
	END IF

	IF (FOUND_MATCH .EQV. .FALSE.) THEN
		WRITE(ERR_MSG,'(A,I4)') ' UNABLE TO FIND MATCHING PERO-FACE FOR iface:', iface
		CALL SHUTDOWN (ERR_MSG)
	END IF
!	WRITE(765,*)IBF, iface, IF2C(iface,3), IF2C(iface,4)
 END DO

 END SUBROUTINE MATCH_PERO_FACES
!----------------------------------------------------------------------------!

