 SUBROUTINE READ_INPUT_FILE

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER, PARAMETER :: finp = 101

 WRITE(*, '(A)', ADVANCE = 'NO') ' READING INPUT FILE... '

 OPEN(finp, FILE = 'SD_input.dat', ACTION = 'READ')

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) PROJECT_NAME
 READ(finp,*)
 READ(finp,*) IS_RESTART

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) VISMODE

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) IS_DEFORMING

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) N

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) INIT_COND

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) ITERMAX, DT, RK_stage

 READ(finp,*)
 READ(finp,*)
 READ(finp,*) TECFREQ, OUT_FORMAT, NRE

 READ(finp,*)
 READ(finp,*)
 READ(finp,*)
 READ(finp,*) rhoinf, uinf, vinf, pinf

 CLOSE(finp)

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''
 
 END SUBROUTINE READ_INPUT_FILE
 !----------------------------------------------------------------------------!


 SUBROUTINE READ_VRT_DATA

 USE setup
 USE MESH2D

 INTEGER, PARAMETER :: fvrt = 101
 INTEGER :: IV, dummy

 CHARACTER (LEN = 50) :: NAMEVRT

 WRITE(*, '(A)', ADVANCE = 'NO') ' READING VRT DATA... '
 
 NAMEVRT = TRIM(PROJECT_NAME)//'.vrt'

 OPEN(fvrt, FILE = TRIM(NAMEVRT), ACTION = 'READ')

 READ(fvrt,*) NVERT

 ALLOCATE(XV(NVERT), YV(NVERT))

 DO IV = 1, NVERT
	READ(fvrt,*) dummy, XV(IV), YV(IV)
 END DO

 CLOSE(fvrt)

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''
 
 END SUBROUTINE READ_VRT_DATA
 !----------------------------------------------------------------------------!


 SUBROUTINE READ_CEL_DATA

 USE setup
 USE MESH2D

 INTEGER, PARAMETER :: fcel = 102
 INTEGER :: IC, dummy

 CHARACTER (LEN = 50) :: NAMECEL

 WRITE(*, '(A)', ADVANCE = 'NO') ' READING CEL DATA... '

 NAMECEL = TRIM(PROJECT_NAME)//'.cel'

 OPEN(fcel, FILE = TRIM(NAMECEL), ACTION = 'READ')

 READ(fcel,*) NCELL

 ALLOCATE(IVCELL(NCELL,4))

 DO IC = 1, NCELL
	READ(fcel,*) dummy, IVCELL(IC,:)
 END DO

 CLOSE(fcel)

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''

 END SUBROUTINE READ_CEL_DATA
!----------------------------------------------------------------------------!


 SUBROUTINE READ_BND_DATA

 USE setup
 USE MESH2D

 INTEGER, PARAMETER :: fbnd = 101
 INTEGER :: IBF, dummy

 CHARACTER (LEN = 50) :: NAMEBND

 WRITE(*, '(A)', ADVANCE = 'NO') ' READING BND DATA... '
 
 NAMEBND = TRIM(PROJECT_NAME)//'.bnd'

 OPEN(fbnd, FILE = TRIM(NAMEBND), ACTION ='READ')

 READ(fbnd,*) NBOUN

 ALLOCATE(IBF2F(NBOUN), IVBOUN(NBOUN,2))

 DO IBF = 1, NBOUN
	READ(fbnd,*) dummy, IBF2F(IBF), IVBOUN(IBF,:), FACEMARKER_BOUN(IBF)
 END DO

 CLOSE(fbnd)

 WRITE(*,'(A)') 'DONE!'
 PRINT*,''

 END SUBROUTINE READ_BND_DATA
!----------------------------------------------------------------------------!
