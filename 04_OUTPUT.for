
!     OUTPUT

!     ===================================================================
      SUBROUTINE OUTPUT(iMAXNOS,iANATYPE,iLOADTYPE,iNODT,rDT,rPI,
     & rMASS,rDAMPCOEFF,rSTIFF,rPZERO,rFREQ,rDAMPRATIO,rWN,rWD,rTN,rTD,
     & rLV,rU,rUD,rUDD,rUDDD,rDURA,iPULSETYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION     
!     ===================================================================
      WRITE(*,*) '04 PROCEDURE: OUTPUT'
      
      OPEN(01,FILE='OUTPUT.TXT')
      
      WRITE(01,*)
      WRITE(01,*) '##_NUMERICAL_ANALYSIS_OUTPUT_##'
      WRITE(01,*)
      WRITE(01,*) '<_SYSTEM_INFORMATION_>'
      WRITE(01,91) rDAMPRATIO
      WRITE(01,92) rWD,rTD
      WRITE(01,*)
      
      WRITE(01,*) '<_LOAD_INFORMATION_>'
      IF(iLOADTYPE==1) CALL OUTPUT_HARMONIC(rPZERO,rFREQ,rPI)
      IF(iLOADTYPE==2) CALL OUTPUT_PULSE(iPULSETYPE,rPZERO,rDURA,rTN)
      IF(iLOADTYPE==3) CALL OUTPUT_GROUND(iMAXNOS)
      IF(iLOADTYPE==4) CALL OUTPUT_GROUND2(iMAXNOS)
      
      WRITE(01,*) '<_ANALYSIS_METHOD_INFORMATION_>'
      IF(iANATYPE==1) WRITE(01,*) 'LINEAR_INTERPOLATION_METHOD'
      IF(iANATYPE==2) WRITE(01,*) 'CENTRAL_DIFFERENCE_METHOD'
      IF(iANATYPE==3) WRITE(01,*) 'NEWMARK`S_METHOD(AVERAGE_ACC.)'
      IF(iANATYPE==4) WRITE(01,*) 'NEWMARK`S_METHOD(LINEAR_ACC.)'
      IF(iANATYPE==5) WRITE(01,*) 'DFT_METHOD'
      IF(iANATYPE==6) WRITE(01,*) 'IMPROVED_DFT_METHOD'
      IF(iANATYPE==7) WRITE(01,*) 'STATE_FORM_METHOD'
      WRITE(01,*)
      
      WRITE(01,*) '<_ANALYSIS_RESULT_>'
      WRITE(01,94)
      iCOUNT=0
      DO 101 i=1,iNODT
        WRITE(01,93) iCOUNT*rDT,rU(i),rUD(i),rUDD(i)
        iCOUNT=iCOUNT+1 
101   CONTINUE        
      
91    FORMAT(' DAMPING_RATIO: ',E16.10)      
92    FORMAT(' NATURAL_FREQUENCY_AND_PERIOD: ',E16.10,3X,E16.10)    
93    FORMAT(3X,E10.4,3X,E16.10,3X,E16.10,3X,E16.10)   
94    FORMAT(3X,'TIME',9X,'DISPLACEMENT',7X,
     & 'VELOCITY',11X,'ACCELERATION')
     
      RETURN
      END SUBROUTINE
!     ===================================================================

!     ===================================================================
      SUBROUTINE OUTPUT_HARMONIC(rPZERO,rFREQ,rPI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================
      WRITE(01,81) rPZERO,rFREQ*2*rPI
      WRITE(01,*)
81    FORMAT(' APPLIED_LOAD: ',E10.4,'*SIN[',E10.4,'*t]') 
      
      RETURN
      END SUBROUTINE     
!     ===================================================================

!     ===================================================================
      SUBROUTINE OUTPUT_PULSE(iPULSETYPE,rPZERO,rDURA,rTN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================
      IF(iPULSETYPE==1) WRITE(01,*) 'RECTANGULAR_PULSE_FORCE'
      IF(iPULSETYPE==2) WRITE(01,*) 'TRIANGULAR_PULSE_FORCE'
      IF(iPULSETYPE==3) WRITE(01,*) 'HALF-SINE_PULSE_FORCE'
      WRITE(01,82) rPZERO,rDURA,rDURA/rTN
      WRITE(01,*)
82    FORMAT(' AMPLITUDE: ',E10.4,3X,'DURATION: ',
     & E10.4,3X,'td/Tn: ',E10.4) 
      
      RETURN
      END SUBROUTINE     
!     ===================================================================

!     ===================================================================
      SUBROUTINE OUTPUT_GROUND(iMAXNOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================
      WRITE(01,*) 'DATA: EL_CENTRO_1940_NORTH_SOUTH_COMPONENT'
      WRITE(01,*)
      
      RETURN
      END SUBROUTINE     
!     ===================================================================

!     ===================================================================
      SUBROUTINE OUTPUT_GROUND2(iMAXNOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================
      WRITE(01,*) 'DATA: EL_CENTRO_1940_NORTH_SOUTH_COMPONENT(2)'
      WRITE(01,*)
      
      RETURN
      END SUBROUTINE     
!     ===================================================================
