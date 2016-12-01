 !---------------------------------------------------------!
 !Program SD_Hyper (Spectral Difference--Hyperbolic)
 !Written independently by Shantanu Bailoor
 !---------------------------------------------------------!

 !---------------------------------------------------------!
 !Code written to obtain high-order accurate simulations of 
 !compressible hyperbolic conservation laws.
 
 !Code based on high-order Spectral Difference (SD) Method
 !Integrates two-dimensional, compressible Euler equations 
 !using high-order Runge-Kutta schemes.
 !---------------------------------------------------------!

 PROGRAM SD_HYPER

 CALL DISP_INIT_MESSAGE

 CALL READ_INPUT_FILE

 CALL READ_VRT_DATA

 CALL READ_CEL_DATA

 CALL READ_BND_DATA

 CALL CONNECTIVITY

 CALL DISP_GEOM_REPORT

 CALL INITSETUP

 CALL ITERATIONS

 CALL FREE_MEMORY
 
 END PROGRAM SD_HYPER
 !----------------------------------------------------------------------------!
