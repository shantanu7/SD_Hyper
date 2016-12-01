 SUBROUTINE WRITE_TEC_DATA_CELLCENTERED

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER, PARAMETER :: ftec = 101
 INTEGER :: IC, IV, is, js, k, nid
 
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: plotQ 
 REAL(KIND=DP) :: xi, eta

 CHARACTER (LEN = 50) :: NAMETEC

 !Create output buffers
 ALLOCATE(plotQ(4,NCELL))
 plotQ(:,:) = 0._DP

 DO IC = 1, NCELL
	xi = 0.5; eta = 0.5 !Coordinates for cell-center of a standard quad.

	DO is = 1, N
		DO js = 1, N
			plotQ(:,IC) = plotQ(:,IC) + Q(:,is,js,IC)*hhval(is,xi)*hhval(js,eta)
		END DO
	END DO
	plotQ(2,IC) = plotQ(2,IC)/plotQ(1,IC)
	plotQ(3,IC) = plotQ(3,IC)/plotQ(1,IC)
	plotQ(4,IC) = (plotQ(4,IC) - 0.5_DP*plotQ(1,IC)* & 
				  (plotQ(2,IC)**2 + plotQ(3,IC)**2))*(gam-1._DP)
 END DO

 WRITE(NAMETEC,'(A,I6.6,A)')'tec2D_SD.',ITER,'.dat'
 OPEN(ftec, NAME=TRIM(NAMETEC), ACTION='WRITE')

 WRITE(ftec,*) 'TITLE = "TECPLOT-COMPATIBLE OUTPUT FILE"'
 WRITE(ftec,*) 'VARIABLES = "X", "Y", "RHO", "U", "V", "P"'
 WRITE(ftec,*) 'ZONE N = ', NVERT, ' E = ', NCELL
 WRITE(ftec,'(A,F)') 'StrandID = 1, SolutionTime = ', ctime
 WRITE(ftec,'(A)') 'DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, VARLOCATION=([3-6]=CELLCENTERED)'

 WRITE(ftec,'(4(2X,E22.13))') XVMP(:), YVMP(:), plotQ(1,:), plotQ(2,:), plotQ(3,:), plotQ(4,:)

 DO IC = 1, NCELL
	WRITE(ftec,'(4(I5,2X))') IVCELL(IC,:)
 END DO

 CLOSE(ftec)

 DEALLOCATE(plotQ)
 !----------------------------------------------------------!


 CONTAINS

 REAL(KIND=DP) FUNCTION hhval(i, xval)

 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), INTENT (IN) :: xval
 REAL(KIND=DP) :: nume, deno

 INTEGER, INTENT (IN) :: i
 INTEGER :: s

 nume = 1._DP; deno = 1._DP

 DO s = 1, N
	IF (s .NE. i) THEN
		nume = nume*(xval-Xs(s))
		deno = deno*(Xs(i)-Xs(s))
	END IF
 END DO 

 hhval = nume/deno
  
 END FUNCTION hhval

 END SUBROUTINE WRITE_TEC_DATA_CELLCENTERED
!----------------------------------------------------------------------------!


 SUBROUTINE WRITE_TEC_DATA_NODAL

 USE MESH2D
 USE setup
 IMPLICIT NONE

 INTEGER, ALLOCATABLE, DIMENSION(:) :: plotQ_fill
 INTEGER, PARAMETER :: ftec = 101
 INTEGER :: IC, IV, is, js, k, nid
 
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: plotQ 
 REAL(KIND=DP) :: xi, eta

 CHARACTER (LEN = 50) :: NAMETEC

 !Create output buffers
 ALLOCATE(plotQ(4,NVERT), plotQ_fill(NVERT))
 plotQ(:,:) = 0._DP; plotQ_fill = 0

 DO IC = 1, NCELL
 	DO k = 1, 4
		IV = IVCELL(IC,k)
		IF (plotQ_fill(IV) .EQ. 0) THEN

			SELECT CASE(k)
				CASE(1)
				xi = 0._DP; eta = 0._DP
				CASE(2)
				xi = 1._DP; eta = 0._DP
				CASE(3)
				xi = 1._DP; eta = 1._DP
				CASE(4)
				xi = 0._DP; eta = 1._DP
			END SELECT

			DO is = 1, N
				DO js = 1, N
					plotQ(:,IV) = plotQ(:,IV) + Q(:,is,js,IC)* &
								  hhval(is,xi)*hhval(js,eta)
				END DO
			END DO

			plotQ(2,IV) = plotQ(2,IV)/plotQ(1,IV)
			plotQ(3,IV) = plotQ(3,IV)/plotQ(1,IV)
			plotQ(4,IV) = (plotQ(4,IV) - 0.5_DP*plotQ(1,IV)* & 
						  (plotQ(2,IV)**2 + plotQ(3,IV)**2))*(gam-1._DP)
			plotQ_fill(IV) = 1

		END IF
	END DO	
 END DO


 WRITE(NAMETEC,'(A,I6.6,A)')'tec2D_SD.',ITER,'.dat'
 OPEN(ftec, NAME=TRIM(NAMETEC), ACTION='WRITE')

 WRITE(ftec,*) 'TITLE = "TECPLOT-COMPATIBLE OUTPUT FILE"'
 WRITE(ftec,*) 'VARIABLES = "X", "Y", "RHO", "U", "V", "P"'
 WRITE(ftec,*) 'ZONE N = ', NVERT, ' E = ', NCELL
 WRITE(ftec,'(A,F)') 'StrandID = 1, SolutionTime = ', ctime
 WRITE(ftec,'(A)') 'DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL'

 WRITE(ftec,'(4(2X,E22.13))') XVMP(:), YVMP(:), plotQ(1,:), plotQ(2,:), plotQ(3,:), plotQ(4,:)

 DO IC = 1, NCELL
	WRITE(ftec,'(4(I5,2X))') IVCELL(IC,:)
 END DO

 CLOSE(ftec)

 DEALLOCATE(plotQ, plotQ_fill)
 !----------------------------------------------------------!


 CONTAINS

 REAL(KIND=DP) FUNCTION hhval(i, xval)

 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), INTENT (IN) :: xval
 REAL(KIND=DP) :: nume, deno

 INTEGER, INTENT (IN) :: i
 INTEGER :: s

 nume = 1._DP; deno = 1._DP

 DO s = 1, N
	IF (s .NE. i) THEN
		nume = nume*(xval-Xs(s))
		deno = deno*(Xs(i)-Xs(s))
	END IF
 END DO 

 hhval = nume/deno
  
 END FUNCTION hhval

 END SUBROUTINE WRITE_TEC_DATA_NODAL
!----------------------------------------------------------------------------!


 SUBROUTINE WRITE_TEC_DATA_NODAL_BINARY

 USE setup
 USE MESH2D
 IMPLICIT NONE

 INTEGER, ALLOCATABLE, DIMENSION(:) :: plotQ_fill
 INTEGER, PARAMETER :: ftec = 101 
 INTEGER :: IC, IV, is, js, k, nid
 
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: plotQ 
 REAL(KIND=DP) :: xi, eta

 CHARACTER (LEN = 50) :: NAMETEC

 !Tecplot binary format related
 !-----------------------------
 INTEGER, PARAMETER :: NVARS = 6, MAX_STRING_LEN = 64 
 INTEGER, ALLOCATABLE, DIMENSION(:) :: itemp
 INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IVCELL_plot

 REAL(KIND=4), PARAMETER :: ZONEMARKER = 299.0, EOHMARKER = 357.0 
 REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: minmax

 CHARACTER (LEN=100) :: string
 
 !Create output buffers
 ALLOCATE(plotQ(4,NVERT), plotQ_fill(NVERT))
 plotQ(:,:) = 0.0; plotQ_fill = 0

 DO IC = 1, NCELL
 	DO k = 1, 4
		IV = IVCELL(IC,k)
		IF (plotQ_fill(IV) .EQ. 0) THEN

			SELECT CASE(k)
				CASE(1)
				xi = 0._DP; eta = 0._DP
				CASE(2)
				xi = 1._DP; eta = 0._DP
				CASE(3)
				xi = 1._DP; eta = 1._DP
				CASE(4)
				xi = 0._DP; eta = 1._DP
			END SELECT

			DO is = 1, N
				DO js = 1, N
					plotQ(:,IV) = plotQ(:,IV) + REAL(Q(:,is,js,IC)* &
								  hhval(is,xi)*hhval(js,eta),KIND=8)
				END DO
			END DO

			plotQ(2,IV) = plotQ(2,IV)/plotQ(1,IV)
			plotQ(3,IV) = plotQ(3,IV)/plotQ(1,IV)
			plotQ(4,IV) = (plotQ(4,IV) - 0.5*plotQ(1,IV)* & 
						  (plotQ(2,IV)**2 + plotQ(3,IV)**2))*(gam-1.0)
			plotQ_fill(IV) = 1

		END IF
	END DO	
 END DO

 !Prepare other variables for file header
 ALLOCATE(minmax(NVARS,2), itemp(MAX_STRING_LEN), IVCELL_plot(NCELL,4))

 !Open output file:
 WRITE(NAMETEC,'(A,I6.6,A)')'tec2D_SD.',ITER,'.plt'		!Binary files will have extension .plt
 OPEN(ftec, NAME=TRIM(NAMETEC), ACCESS='STREAM', FORM='UNFORMATTED', ACTION='WRITE')

 WRITE(ftec) '#!TDV112'	!Magic number (Tecplot Version)

 !Title of file
 itemp(1) = 1			!Integer value of 1
 itemp(2) = 0			!File type (0=full, 1=grid, 2=solution)
 WRITE(ftec) itemp(:2)

 !Start writing file header
 string = 'SD_Hyper Data' 	!Title of data set
 CALL WRITE_TEC_STRING(ftec, string)

 !NVARS and variable names
 WRITE(ftec) NVARS

 string = "X"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "Y"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "RHO"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "U"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "V"
 CALL WRITE_TEC_STRING(ftec, string)

 string = "P"
 CALL WRITE_TEC_STRING(ftec, string)

 WRITE(ftec) ZONEMARKER

 string = 'ZONE 001'
 CALL WRITE_TEC_STRING(ftec, string)

 itemp(1) =-1		!Parent Zone
 itemp(2) = 0		!Strand ID
 WRITE(ftec) itemp(:2)

 WRITE(ftec) REAL(ctime,KIND=8)

 !Dataset format information
 itemp(1)  = -1		!Zone color
 itemp(2)  =  3		!Zone type (3=FEQUADRILATERAL)
 itemp(3)  =  0		!Data packing (0=block, 1=point)
 itemp(4)  =  0		!Var location (0=nodal, 1=specify)
 itemp(5)  =  0		!0=nodal, 1=cellcentered
 itemp(6)  = NVERT	!Number of vertices
 itemp(7)  = NCELL	!Number of cells
 itemp(8)  =  0		!Icell dim (for i-, j-, k- type data)
 itemp(9)  =  0		!Jcell dim (for i-, j-, k- type data)
 itemp(10) =  0		!Kcell dim (for i-, j-, k- type data)
 itemp(11) =  0		!Other auxiliary data

 WRITE(ftec) itemp(:11)

 WRITE(ftec) EOHMARKER

 !Write Zones
 WRITE(ftec) ZONEMARKER

 itemp(1:NVARS) = 2	!Variable format: 1=float, 2=double...
 itemp(NVARS+1) = 0	!Has passive variables?
 itemp(NVARS+2) = 0	!Has variable sharing?
 itemp(NVARS+3) =-1	!Zone to share connectivity (-1=no)

 WRITE(ftec) itemp(:NVARS+3)

 minmax(1,1) = MINVAL(XVMP); 		minmax(1,2) = MAXVAL(XVMP);
 minmax(2,1) = MINVAL(YVMP); 		minmax(2,2) = MAXVAL(YVMP);
 minmax(3,1) = MINVAL(plotQ(1,:)); 	minmax(3,2) = MAXVAL(plotQ(1,:));
 minmax(4,1) = MINVAL(plotQ(2,:)); 	minmax(4,2) = MAXVAL(plotQ(2,:));
 minmax(5,1) = MINVAL(plotQ(3,:)); 	minmax(5,2) = MAXVAL(plotQ(3,:));
 minmax(6,1) = MINVAL(plotQ(4,:)); 	minmax(6,2) = MAXVAL(plotQ(4,:));

 DO k = 1, NVARS
 	WRITE(ftec) REAL(minmax(k,1),KIND=8)
 	WRITE(ftec) REAL(minmax(k,2),KIND=8)
 END DO

 !Start writing variables
 WRITE(ftec) REAL(XVMP(:NVERT),KIND=8)
 WRITE(ftec) REAL(YVMP(:NVERT),KIND=8)
 WRITE(ftec) REAL(plotQ(1,:NVERT),KIND=8)
 WRITE(ftec) REAL(plotQ(2,:NVERT),KIND=8)
 WRITE(ftec) REAL(plotQ(3,:NVERT),KIND=8)
 WRITE(ftec) REAL(plotQ(4,:NVERT),KIND=8)

 !Write connectivity
 DO IC = 1, NCELL
	DO k = 1, 4
		IVCELL_plot(IC,k) = IVCELL(IC,k) - 1
	END DO
 END DO

 DO IC = 1, NCELL
	 WRITE(ftec) IVCELL_plot(IC,:)
 END DO

 CLOSE(ftec)

 DEALLOCATE(plotQ, plotQ_fill)
 DEALLOCATE(minmax, itemp, IVCELL_plot)
 
 !----------------------------------------------------------!


 CONTAINS

 REAL(KIND=DP) FUNCTION hhval(i, xval)

 USE setup
 IMPLICIT NONE

 REAL(KIND=DP), INTENT (IN) :: xval
 REAL(KIND=DP) :: nume, deno

 INTEGER, INTENT (IN) :: i
 INTEGER :: s

 nume = 1._DP; deno = 1._DP

 DO s = 1, N
	IF (s .NE. i) THEN
		nume = nume*(xval-Xs(s))
		deno = deno*(Xs(i)-Xs(s))
	END IF
 END DO 

 hhval = nume/deno
  
 END FUNCTION hhval
 !----------------------------------------------------------!

 SUBROUTINE WRITE_TEC_STRING(ftec, string)

 IMPLICIT NONE

 INTEGER, INTENT(IN) :: ftec
 INTEGER, ALLOCATABLE, DIMENSION (:) :: int_str
 INTEGER :: i, strlen

 CHARACTER(LEN=100), INTENT(IN) :: string

 strlen = LEN_TRIM(string)

 ALLOCATE(int_str(strlen+1))

 DO i = 1, strlen
	int_str(i) = ICHAR(string(i:i))
 END DO
 int_str(strlen+1) = 0

 WRITE(ftec) int_str(1:strlen+1)

 DEALLOCATE(int_str)
 
 END SUBROUTINE WRITE_TEC_STRING

 END SUBROUTINE WRITE_TEC_DATA_NODAL_BINARY
!----------------------------------------------------------------------------!
