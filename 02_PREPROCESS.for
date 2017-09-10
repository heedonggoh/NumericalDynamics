
!     SUBROUTINE PREPROCESS

!     ===================================================================      
      SUBROUTINE PREPROCESS(iMAXNOS,rPI,rMASS,rDAMPCOEFF,rSTIFF,rDT,
     & iNODT,iLOADTYPE,rPZERO,rFREQ,rUDDG,
     & rLV,rWN,rDAMPRATIO,rWD,rTN,rTD,rDURA,iPULSETYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rUDDG(iMAXNOS)
!     rLV = LOAD VECTOR  
!     rUDDG = GROUND MOTION
!     ===================================================================      
      WRITE(*,*) '02 PROCEDURE: PREPROCESS'

!     SAVE LOAD VECTOR      
      IF(iLOADTYPE==1) CALL PREPROCESS_HARMONIC(iMAXNOS,rPI,rDT,iNODT,
     & rPZERO,rFREQ,rLV)
      IF(iLOADTYPE==2) CALL PREPROCESS_PULSE(iMAXNOS,rPI,rDT,iNODT,
     & rPZERO,rDURA,rLV,iPULSETYPE)
      IF(iLOADTYPE==3) CALL PREPROCESS_GROUND(iMAXNOS,rDT,iNODT,rLV,
     & rUDDG,rMASS)
      IF(iLOADTYPE==4) CALL PREPROCESS_GROUND(iMAXNOS,rDT,iNODT,rLV,
     & rUDDG,rMASS)
     
!     SAVE NATURAL AND DYNAMIC ANGULAR FREQUENCY AND DAMPING RATIO      
      rWN=DSQRT(rSTIFF/rMASS)
      rDAMPRATIO=rDAMPCOEFF/(2.0D0*DSQRT(rSTIFF*rMASS))
      rWD=rWN*DSQRT(1.0D0-rDAMPRATIO**2)
      rTN=2*rPI/rWN
      rTD=2*rPI/rWD
      
!     RESONANCE CHECK
      IF(rDAMPCOEFF==0.0D0.AND.rWN==rFREQ) THEN
        WRITE(*,*) 'WARNING: RESONANCE EXPECTED'
        STOP
      END IF        
      
      RETURN
      END SUBROUTINE
!     ===================================================================      

!     ===================================================================      
      SUBROUTINE PREPROCESS_HARMONIC(iMAXNOS,rPI,rDT,iNODT,rPZERO,
     & rFREQ,rLV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS)
!     rLV = LOAD VECTOR  
!     ===================================================================      
      WRITE(*,*) '   LOAD_TYPE: HARMONIC LOADING'
      
      iCOUNT=0
      rLV(:)=0.0D0
      DO 101 i=1,iNODT
        rLV(i)=rPZERO*DSIN(2.0D0*rPI*rFREQ*rDT*iCOUNT)
        iCOUNT=iCOUNT+1
101   CONTINUE      

      RETURN
      END SUBROUTINE
!     ===================================================================            

!     ===================================================================      
      SUBROUTINE PREPROCESS_PULSE(iMAXNOS,rPI,rDT,iNODT,
     & rPZERO,rDURA,rLV,iPULSETYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS)
!     rLV = LOAD VECTOR  
!     ===================================================================      
      WRITE(*,*) '   LOAD_TYPE: PULSE LOADING'
      
      iDURA=NINT(rDURA/rDT)
      rLV(:)=0.0D0
      IF(iPULSETYPE==1) THEN
        DO 101 i=1,iDURA
            rLV(i)=rPZERO
101     CONTINUE
      ELSE IF(iPULSETYPE==2) THEN
        iCOUNT=0
        rA=2.0D0*rPZERO/rDURA
        DO 201 i=1,iDURA/2
            rLV(i)=rA*iCOUNT*rDT
            rLV(i+iDURA/2)=-rA*iCOUNT*rDT+rPZERO
            iCOUNT=iCOUNT+1
201     CONTINUE  
      ELSE IF(iPULSETYPE==3) THEN
        iCOUNT=0
        DO 301 i=1,iDURA
            rLV(i)=rPZERO*DSIN(rPI*iCOUNT*rDT/rDURA)
            iCOUNT=iCOUNT+1
301     CONTINUE          
      END IF            
      
      RETURN
      END SUBROUTINE
!     ===================================================================     

!     ===================================================================      
      SUBROUTINE PREPROCESS_GROUND(iMAXNOS,rDT,iNODT,rLV,
     & rUDDG,rMASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rUDDG(iMAXNOS)
!     rLV = LOAD VECTOR  
!     ===================================================================      
      WRITE(*,*) '   LOAD_TYPE: EARTHQUAKE LOADING'
  
      rLV(:)=0.0D0
      DO 101 i=1,iNODT
        rLV(i)=-rMASS*rUDDG(i)*9.81D0
101   CONTINUE      

      RETURN
      END SUBROUTINE
!     ===================================================================            
       
