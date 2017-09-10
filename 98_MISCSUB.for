
!     MISC. SUBROUTINES (CREATED 2010.11.17)
	
***** Matrix * Vector *****
	SUBROUTINE MATMULTI(rMATRIX,iSIZE,rINPUT,iRSIZE,iCSIZE)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
	DIMENSION rMATRIX(iSIZE,iSIZE),rINPUT(iSIZE),rOUTPUT(iSIZE)
	rOUTPUT(:)=0.0D0
		
	DO 2001 i=1,iRSIZE
	  DO 2002 j=1,iCSIZE
	    rOUTPUT(i)=rOUTPUT(i)+rMATRIX(i,j)*rINPUT(j)
2002    CONTINUE
2001  CONTINUE	    
      
      rINPUT(:)=0.0D0
      DO 2003 i=1,iSIZE
        rINPUT(i)=rOUTPUT(i)
2003  CONTINUE	           

	RETURN
	END SUBROUTINE
	
***** Matrix * Matrix *****
      SUBROUTINE MATMULT(rA,iAHEIGHT,iAWIDTH,rB,iBHEIGHT,iBWIDTH,rP)	
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION rA(iAHEIGHT,iAWIDTH),rB(iBHEIGHT,iBWIDTH),
     & rP(iAHEIGHT,iBWIDTH)
      rP(:,:)=0.0D0
      
      DO 3001 i=1,iAHEIGHT
        DO 3002 j=1,iBWIDTH
            DO 3003 k=1,iAWIDTH
                rP(i,j)=rP(i,j)+rA(i,k)*rB(k,j)
3003        CONTINUE
3002    CONTINUE

3001  CONTINUE    

      RETURN
      END SUBROUTINE   

***** COMPLEX Matrix * COMPLEX Matrix *****
      SUBROUTINE CMATMULT(rA,iAHEIGHT,iAWIDTH,rB,iBHEIGHT,iBWIDTH,rP)	
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DOUBLE COMPLEX rA(iAHEIGHT,iAWIDTH),rB(iBHEIGHT,iBWIDTH),
     & rP(iAHEIGHT,iBWIDTH)
      rP(:,:)=DCMPLX(0.0D0,0.0D0)
      
      DO 3001 i=1,iAHEIGHT
        DO 3002 j=1,iBWIDTH
            DO 3003 k=1,iAWIDTH
                rP(i,j)=rP(i,j)+rA(i,k)*rB(k,j)
3003        CONTINUE
3002    CONTINUE
3001  CONTINUE    

      RETURN
      END SUBROUTINE   



***** Transpose *****
      SUBROUTINE TRANS(rA,iAHEIGHT,iAWIDTH,rP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION rA(iAHEIGHT,iAWIDTH),rP(iAWIDTH,iAHEIGHT)
      rP(:,:)=0.0D0
      
      DO 101 i=1,iAHEIGHT
        DO 102 j=1,iAWIDTH
            rP(j,i)=rA(i,j)
102     CONTINUE
101   CONTINUE            
      
      RETURN
      END SUBROUTINE
      
***** Additional Matrix Multiplication [A^T]*[B]*[A] *****
      SUBROUTINE MULTATBA(rA,iAHEIGHT,iAWIDTH,rB,rP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION rA(iAHEIGHT,iAWIDTH),rB(iAHEIGHT,iAHEIGHT),
     & rAT(iAWIDTH,iAHEIGHT),rP(iAWIDTH,iAWIDTH),rATB(iAWIDTH,iAHEIGHT)
      rP(:,:)=0.0D0
      rAT(:,:)=0.0D0
      
      CALL TRANS(rA,iAHEIGHT,iAWIDTH,rAT)
      CALL MATMULT(rAT,iAWIDTH,iAHEIGHT,rB,iAHEIGHT,iAHEIGHT,rATB)
      CALL MATMULT(rATB,iAWIDTH,iAHEIGHT,rA,iAHEIGHT,iAWIDTH,rP)
      
      RETURN
      END SUBROUTINE          
      
***** Additional Matrix Multiplication [A]*[B]*[A^T] *****
      SUBROUTINE MULTABAT(rA,iAHEIGHT,iAWIDTH,rB,rP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION rA(iAHEIGHT,iAWIDTH),rB(iAWIDTH,iAWIDTH),
     & rAT(iAWIDTH,iAHEIGHT),rP(iAHEIGHT,iAHEIGHT),rAB(iAHEIGHT,iAWIDTH)
      rP(:,:)=0.0D0
      rAT(:,:)=0.0D0
      
      CALL TRANS(rA,iAHEIGHT,iAWIDTH,rAT)
      CALL MATMULT(rA,iAHEIGHT,iAWIDTH,rB,iAWIDTH,iAWIDTH,rAB)
      CALL MATMULT(rAB,iAHEIGHT,iAWIDTH,rAT,iAWIDTH,iAHEIGHT,rP)
      
      RETURN
      END SUBROUTINE          
      
      
      
***** ALLOCATE DIMENSION *****
!     1st ORDER DOUBLE PRECISION     
      SUBROUTINE rALLOC1ST(rA,iHEIGHT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      REAL*8, DIMENSION(:), ALLOCATABLE :: rA
      
      ALLOCATE(rA(iHEIGHT), STAT=iALLOCATESTATUS)
      IF (iALLOCATESTATUS.NE.0) THEN
        WRITE(*,*) "ALLOCATE ERROR: OVERFLOW ON ARRAY SIZE CALCULATION"
        STOP
      END IF
      rA(:)=0.0D0
      
      RETURN
      END SUBROUTINE
!    2nd ORDER DOUBLE PRECISION    
      SUBROUTINE rALLOC2ND(rA,iHEIGHT,iWIDTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: rA
      
      ALLOCATE(rA(iHEIGHT,iWIDTH), STAT=iALLOCATESTATUS)
      IF (iALLOCATESTATUS.NE.0) THEN
        WRITE(*,*) "ALLOCATE ERROR: OVERFLOW ON ARRAY SIZE CALCULATION"
        STOP
      END IF
      rA(:,:)=0.0D0
      
      RETURN
      END SUBROUTINE
!    3rd ORDER DOUBLE PRECISION      
      SUBROUTINE rALLOC3RD(rA,iHEIGHT,iWIDTH,iDEPTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: rA
      
      ALLOCATE(rA(iHEIGHT,iWIDTH,iDEPTH), STAT=iALLOCATESTATUS)
      IF (iALLOCATESTATUS.NE.0) THEN
        WRITE(*,*) "ALLOCATE ERROR: OVERFLOW ON ARRAY SIZE CALCULATION"
        STOP
      END IF
      rA(:,:,:)=0.0D0
      
      RETURN
      END SUBROUTINE     
!     1st ORDER INTEGER     
      SUBROUTINE iALLOC1ST(iA,iHEIGHT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      INTEGER, DIMENSION(:), ALLOCATABLE :: iA
      
      ALLOCATE(iA(iHEIGHT), STAT=iALLOCATESTATUS)
      IF (iALLOCATESTATUS.NE.0) THEN
        WRITE(*,*) "ALLOCATE ERROR: OVERFLOW ON ARRAY SIZE CALCULATION"
        STOP
      END IF
      iA(:)=0
      
      RETURN
      END SUBROUTINE
!    2nd ORDER INTEGER    
      SUBROUTINE iALLOC2ND(iA,iHEIGHT,iWIDTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iA
      
      ALLOCATE(iA(iHEIGHT,iWIDTH), STAT=iALLOCATESTATUS)
      IF (iALLOCATESTATUS.NE.0) THEN
        WRITE(*,*) "ALLOCATE ERROR: OVERFLOW ON ARRAY SIZE CALCULATION"
        STOP
      END IF
      iA(:,:)=0
      
      RETURN
      END SUBROUTINE
!    3rd ORDER INTEGER      
      SUBROUTINE iALLOC3RD(iA,iHEIGHT,iWIDTH,iDEPTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: iA
      
      ALLOCATE(iA(iHEIGHT,iWIDTH,iDEPTH), STAT=iALLOCATESTATUS)
      IF (iALLOCATESTATUS.NE.0) THEN
        WRITE(*,*) "ALLOCATE ERROR: OVERFLOW ON ARRAY SIZE CALCULATION"
        STOP
      END IF
      iA(:,:,:)=0
      
      RETURN
      END SUBROUTINE        
      
***** PRINT VECTOR *****    
      SUBROUTINE DPRINTVECTOR(rV,iSIZE,iPRINTSIZE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DIMENSION rV(iSIZE)
      
      OPEN(01,FILE='PRINTVECTOR.TXT')
      
      WRITE(01,99) iPRINTSIZE
      DO 101 i=1,iPRINTSIZE
        WRITE(01,98) rV(i)
101   CONTINUE        
     
99    FORMAT('VECTOR_SIZE: ',I6)
98    FORMAT(E16.10)     
      
      RETURN
      END SUBROUTINE        
                   
***** PRINT MATRIX *****
      SUBROUTINE DPRINTMAT(rMAT,iHEIGHT,iWIDTH,iPH,iPW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION rMAT(iHEIGHT,iWIDTH)
      
      OPEN(77,FILE='MATRIX.txt')
      
      WRITE(77,*)
      WRITE(77,"('MATRIX_SIZE:',I4,' *',I4)") iPH,iPW
      WRITE(77,*)
      
      DO 101 i=1,iPH
        DO 102 j=1,iPW
            WRITE(77,"(E16.10,3X)",ADVANCE="NO") rMAT(i,j)
102     CONTINUE
        WRITE(77,*)
101   CONTINUE
      
      RETURN
      END SUBROUTINE             
      
      
***** PRINT COMPLEX VECTOR *****    
      SUBROUTINE CPRINTVECTOR(rV,iSIZE,iPRINTSIZE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)	
      DOUBLE COMPLEX rV(iSIZE)
      
      OPEN(01,FILE='COMPLEXVECTOR.TXT')
      
      WRITE(01,99) iPRINTSIZE
      DO 101 i=1,iPRINTSIZE
        WRITE(01,98) DREAL(rV(i)),DIMAG(rV(i))
101   CONTINUE        
     
99    FORMAT('COMPLEX_VECTOR_SIZE: ',I6)
98    FORMAT(E16.10,1X,E16.10)     
      
      RETURN
      END SUBROUTINE        
                   
***** PRINT COMPLEX MATRIX *****
      SUBROUTINE CPRINTMAT(rMAT,iHEIGHT,iWIDTH,iPH,iPW)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE COMPLEX rMAT(iHEIGHT,iWIDTH)
      
      OPEN(77,FILE='COMPLEXMATRIX.txt')
      
      WRITE(77,*)
      WRITE(77,"('COMPLEX_MATRIX_SIZE:',I4,' *',I4)") iPH,iPW
      WRITE(77,*)
      
      DO 101 i=1,iPH
        DO 102 j=1,iPW
            WRITE(77,"(E16.10,'+i(',E16.10,')',3X)",ADVANCE="NO") 
     &       DREAL(rMAT(i,j)),DIMAG(rMAT(i,j))
102     CONTINUE
        WRITE(77,*)
101   CONTINUE
      
      RETURN
      END SUBROUTINE                   
