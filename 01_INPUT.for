
!     SUBROUTINE INPUT

!     ===================================================================
      SUBROUTINE INPUT(iMAXNOS,iANATYPE,rMASS,rDAMPCOEFF,rSTIFF,rDT,
     & iNODT,iLOADTYPE,rPZERO,rFREQ,rU,rUD,iPULSETYPE,rDURA,rUDDG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rU(iMAXNOS),rUD(iMAXNOS),rUDDG(iMAXNOS)
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDDG = GROUND MOTION
      
!     ===================================================================
      WRITE(*,*) '01 PROCEDURE: INPUT'
      
      OPEN(01,FILE='INPUT.TXT')

!     READ 
!       ANALYSIS_TYPE
!       SIZE_OF_DELTA_t, NUMBER_OF_DELTA_t 
!       MASS, DAMPING_COEFFICIENT, STIFFNESS
!       INITIAL CONDITIONS
!       LOAD_TYPE      
      READ(01,*)
      READ(01,*)
      READ(01,*) iANATYPE
      READ(01,*)
      READ(01,*) rDT,iNODT
      READ(01,*)
      READ(01,*) rMASS,rDAMPCOEFF,rSTIFF
      READ(01,*)
      READ(01,*) rU(1),rUD(1)
      READ(01,*)
      READ(01,*) iLOADTYPE
 
 !    READ LOAD INFORMATION     
      IF(iLOADTYPE==1) CALL INPUT_HARMONIC(rPZERO,rFREQ)
      IF(iLOADTYPE==2) CALL INPUT_PULSE(iPULSETYPE,rPZERO,rDURA)
      IF(iLOADTYPE==3) CALL INPUT_GROUND(iMAXNOS,rUDDG)
      IF(iLOADTYPE==4) CALL INPUT_GROUND2(iMAXNOS,rUDDG)
      
      RETURN
      END SUBROUTINE
!     ===================================================================      

!     ===================================================================      
      SUBROUTINE INPUT_HARMONIC(rPZERO,rFREQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================      
      WRITE(*,*) '   LOAD TYPE: HARMONIC LOADING'
!     READ
!       MAGNITUDE OF THE APPLIED FORCE
!       FREQUENCY OF THE APPLIED FORCE      
      READ(01,*)
      READ(01,*) rPZERO,rFREQ
      
      RETURN
      END SUBROUTINE
!     ===================================================================      

!     ===================================================================      
      SUBROUTINE INPUT_PULSE(iPULSETYPE,rPZERO,rDURA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================      
      WRITE(*,*) '   LOAD TYPE: PULSE LOADING'
!     READ
!       PULSE TYPE (1=RECTANGULAR, 2=TRIANGULAR, 3=HALF-SINE)
!       MAGNITUDE OF THE APPLIED FORCE
!       DURATION TIME OF THE APPLIED FORCE      
      READ(01,*)
      READ(01,*) iPULSETYPE,rPZERO,rDURA
      
      RETURN
      END SUBROUTINE
!     ===================================================================      

!     ===================================================================      
      SUBROUTINE INPUT_GROUND(iMAXNOS,rUDDG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rUDDG(iMAXNOS)
!     ===================================================================      
      WRITE(*,*) '   LOAD TYPE: EARTHQUAKE LOADING'
!     READ
!       GROUND MOTION
!     < ANALYSIS RESOLUTION > (SIZE_OF_Delta_t, NUMBER_OF_Delta_ts)
!                              0.02	1560
      
      OPEN(02,FILE='ELCENTRO.DAT')
      READ(02,*)
      READ(02,*)
      READ(02,*)
      READ(02,*)
      READ(02,*)
      READ(02,*)
      DO 101 i=1,195
        READ(02,*) rA1,rA2,rA3,rA4,rA5,rA6,rA7,rA8
        iCOUNT=(i-1)*8
        rUDDG(iCOUNT+1)=rA1
        rUDDG(iCOUNT+2)=rA2
        rUDDG(iCOUNT+3)=rA3
        rUDDG(iCOUNT+4)=rA4
        rUDDG(iCOUNT+5)=rA5
        rUDDG(iCOUNT+6)=rA6
        rUDDG(iCOUNT+7)=rA7
        rUDDG(iCOUNT+8)=rA8
101   CONTINUE 

      CLOSE(02)
      RETURN
      END SUBROUTINE
!     ===================================================================      

!     ===================================================================      
      SUBROUTINE INPUT_GROUND2(iMAXNOS,rUDDG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rUDDG(iMAXNOS)
      DIMENSION TEMP(iMAXNOS)
!     ===================================================================      
      WRITE(*,*) '   LOAD TYPE: EARTHQUAKE LOADING'
!     READ
!       GROUND MOTION
      
      OPEN(02,FILE='ELCENTRO.DAT')
      READ(02,*)
      READ(02,*)
      READ(02,*)
      READ(02,*)
      READ(02,*)
      READ(02,*)
      rUDDG(:)=0.0D0
      TEMP(:)=0.0D0
      DO 101 i=1,195
        READ(02,*) rA1,rA2,rA3,rA4,rA5,rA6,rA7,rA8
        iCOUNT=(i-1)*8
        TEMP(iCOUNT+1)=rA1
        TEMP(iCOUNT+2)=rA2
        TEMP(iCOUNT+3)=rA3
        TEMP(iCOUNT+4)=rA4
        TEMP(iCOUNT+5)=rA5
        TEMP(iCOUNT+6)=rA6
        TEMP(iCOUNT+7)=rA7
        TEMP(iCOUNT+8)=rA8
101   CONTINUE

      DO i=1,1560
        rUDDG(2*i-1)=TEMP(i)
        rUDDG(2*i)=(TEMP(i)+TEMP(i+1))*0.5D0
      END DO

      CLOSE(02)
      RETURN
      END SUBROUTINE
!     ===================================================================      
