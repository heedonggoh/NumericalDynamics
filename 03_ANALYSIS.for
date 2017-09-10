
!     SUBROUTINE ANALYSIS

!     ===================================================================
      SUBROUTINE ANALYSIS(iMAXNOS,rSTIFF,rDT,iNODT,iLOADTYPE,iANATYPE,
     & rLV,rWN,rDAMPRATIO,rWD,rU,rUD,rUDD,rMASS,rDAMPCOEFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION     
!     ===================================================================
      WRITE(*,*) '03 PROCEDURE: ANALYSIS'
      
      IF(iANATYPE==1) 
     & CALL ANALYSIS_LINEAR(iMAXNOS,rSTIFF,rDT,iNODT,rLV,rWN,rDAMPRATIO,
     &  rWD,rU,rUD,rUDD,rMASS)
      IF(iANATYPE==2) 
     & CALL ANALYSIS_CDM(iMAXNOS,rSTIFF,rDT,iNODT,rLV,rWN,rU,rUD,rUDD,
     &  rMASS,rDAMPCOEFF)
      IF(iANATYPE==3.OR.IANATYPE==4) 
     & CALL ANALYSIS_NEWMARK(iMAXNOS,rSTIFF,rDT,iNODT,rLV,rWN,rU,rUD,
     &  rUDD,rMASS,rDAMPCOEFF,iANATYPE)
      IF(iANATYPE==5) 
     & CALL ANALYSIS_DFT(iMAXNOS,rSTIFF,rDT,iNODT,rLV,rWN,rWD,rU,rUD,
     &  rUDD,rMASS,rDAMPCOEFF,rDAMPRATIO,iANATYPE)
      IF(iANATYPE==6) 
     & CALL ANALYSIS_IMPROVEDDFT(iMAXNOS,rSTIFF,rDT,iNODT,rLV,rWN,rWD,
     &  rU,rUD,rUDD,rMASS,rDAMPCOEFF,rDAMPRATIO,iANATYPE)
      IF(iANATYPE==7) 
     & CALL ANALYSIS_STATEFORM(iMAXNOS,rSTIFF,rDT,iNODT,rLV,rWN,rWD,rU,
     &  rUD,rUDD,rMASS,rDAMPCOEFF,rDAMPRATIO,iANATYPE)
     
      RETURN
      END SUBROUTINE
!     ===================================================================

!     ===================================================================
      SUBROUTINE ANALYSIS_LINEAR(iMAXNOS,rSTIFF,rDT,iNODT,
     & rLV,rWN,rDAMPRATIO,rWD,rU,rUD,rUDD,rMASS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION     
!     ===================================================================
      WRITE(*,*) '   ANALYSIS TYPE: LINEAR_INTERPOLATION'

!     COEFFICIENTS IN RECURRENCE FORMULAS
      rDAMP=rDAMPRATIO
      rEXP=DEXP(-rDAMP*rWN*rDT)
      rSIN=DSIN(rWD*rDT)
      rCOS=DCOS(rWD*rDT)
      rROOT=DSQRT(1-rDAMP**2)
      rWNDT=rWN*rDT

      rA=rEXP*(rDAMP*rSIN/rROOT+rCOS)
      rB=rEXP*rSIN/rWD
      rC=(1.0D0/rSTIFF)
     & *(2*rDAMP/rWNDT+rEXP*(((1.0D0-2.0D0*rDAMP**2)/rWNDT
     & -rDAMP/rROOT)*rSIN-(1.0D0+2.0D0*rDAMP/rWNDT)*rCOS))
      rD=(1.0D0/rSTIFF)
     & *(1.0D0-2.0D0*rDAMP/rWNDT
     & +rEXP*(((2.0D0*rDAMP**2-1.0D0)/rWNDT)*rSIN
     & +(2.0D0*rDAMP/rWNDT)*rCOS))
      
      rAD=-rEXP*(rWN*rSIN/rROOT)
      rBD=rEXP*(rCOS-rDAMP*rSIN/rROOT)
      rCD=(1.0D0/rSTIFF)
     & *(-1.0D0/rDT+rEXP*((rWN/rROOT+rDAMP/(rDT*rROOT))*rSIN
     & +rCOS/rDT))
      rDD=(1.0D0/(rSTIFF*rDT))*(1.0D0-rEXP*(rDAMP*RSIN/rROOT+rCOS))
      
!     DISPLACEMENT AND VELOCITY
      DO 101 i=1,iNODT-1
        rU(i+1)=rA*rU(i)+rB*rUD(i)+rC*rLV(i)+rD*rLV(i+1)
        rUD(i+1)=rAD*rU(i)+rBD*rUD(i)+rCD*rLV(i)+rDD*rLV(i+1)
101   CONTINUE

!     ACCELERATION
      DO 102 i=1,iNODT
        rUDD(i)=(-rDAMP*rUD(i)-rSTIFF*rU(i)+rLV(i))/rMASS
102   CONTINUE        
      
      RETURN
      END SUBROUTINE
!     ===================================================================     

!     ===================================================================
      SUBROUTINE ANALYSIS_CDM(iMAXNOS,rSTIFF,rDT,iNODT,
     & rLV,rWN,rU,rUD,rUDD,rMASS,rDAMPCOEFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION
!     ===================================================================
      WRITE(*,*) '   ANALYSIS TYPE: CENTRAL_DIFFERENCE_METHOD'

!     BLOW UP CHECK
      IF(rDT*rWN>=2.0D0) THEN
        WRITE(*,*) 'CDM STABILITY CONDITION VIOLATION'
        STOP
      END IF
      
!     INITIAL CALCULATIONS
      rDAMP=rDAMPCOEFF
      rDT2=rDT**2
      rUDD(1)=(rLV(1)-rDAMP*rUD(1)-rSTIFF*rU(1))/rMASS
      rUA=rU(1)-rDT*rUD(1)-rDT2*rUDD(1)/2.0D0
      rKHAT=rMASS/rDT2+rDAMP/(2.0D0*rDT)
      rA=rMASS/rDT2-rDAMP/(2.0D0*rDT)
      rB=rSTIFF-2.0D0*rMASS/rDT2
      
!     CALCULATIONS FOR TIME STEP 1
      rPHAT=rLV(1)-rA*rUA-rB*rU(1)
      rU(2)=rPHAT/rKHAT
      
!     REPETITION FOR THE NEXT TIME STEP
      DO 101 i=2,iNODT-1
        rPHAT=rLV(i)-rA*rU(i-1)-rB*rU(i)
        rU(i+1)=rPHAT/rKHAT
        rUD(i)=(rU(i+1)-rU(i-1))/(2.0D0*rDT)
        rUDD(i)=(rU(i+1)-2.0D0*rU(i)+rU(i-1))/rDT2  
101   CONTINUE          
      
      RETURN
      END SUBROUTINE
!     ===================================================================       

!     ===================================================================
      SUBROUTINE ANALYSIS_NEWMARK(iMAXNOS,rSTIFF,rDT,iNODT,
     & rLV,rWN,rU,rUD,rUDD,rMASS,rDAMPCOEFF,iANATYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION
!     ===================================================================
      WRITE(*,*) '   ANALYSIS TYPE: NEWMARK`S_METHOD'

      IF(iANATYPE==3) THEN
        rGAM=0.5D0
        rBET=0.25
        WRITE(*,*) '                  (AVERAGE ACCELERATION METHOD)'
      ELSE IF(iANATYPE==4) THEN
        rGAM=0.5D0
        rBET=1.0D0/6.0D0
        WRITE(*,*) '                  (LINEAR ACCELERATION METHOD)'
      END IF  

!     BLOW UP CHECK
      IF(rDT*rWN*DSQRT(2*rGAM-4*rBET)>2.0D0) THEN
        WRITE(*,*) 'NEWMARK`S METHOD STABILITY CONDITION VIOLATION'
        STOP
      END IF
      
!     INITIAL CALCULATIONS
      rBETDT=rBET*rDT
      rBETDT2=rBET*(rDT**2)
      rDAMP=rDAMPCOEFF
      
      rUDD(1)=(rLV(1)-rDAMP*rUD(1)-rSTIFF*rU(1))/rMASS      
      rKHAT=rSTIFF+rGAM*rDAMP/rBETDT+rMASS/rBETDT2
      rA=rMASS/rBETDT+rGAM*rDAMP/rBET
      rB=rMASS/(2.0D0*rBET)+rDT*rDAMP*(rGAM/(2.0D0*rBET)-1.0D0)
      
!     CALCULATIONS FOR EACH TIME STEP
      DO 101 i=1,iNODT-1
        rDPHAT=rLV(i+1)-rLV(i)+rA*rUD(i)+rB*rUDD(i)      
        rDU=rDPHAT/rKHAT
        rDUD=rGAM*rDU/(rBET*rDT)
     &   -rGAM*rUD(i)/rBET+rDT*rUDD(i)*(1.0D0-rGAM/(2.0D0*rBET))
        rDUDD=rDU/rBETDT2-rUD(i)/rBETDT-rUDD(i)/(2.0D0*rBET)
        rU(i+1)=rU(i)+rDU
        rUD(i+1)=rUD(i)+rDUD
        rUDD(i+1)=rUDD(i)+rDUDD
101   CONTINUE        
      
      RETURN
      END SUBROUTINE
!     ===================================================================
      
!     ===================================================================
      SUBROUTINE ANALYSIS_DFT(iMAXNOS,rSTIFF,rDT,iNODT,
     & rLV,rWN,rWD,rU,rUD,rUDD,rMASS,rDAMPCOEFF,rDAMPRATIO,iANATYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION
!     * LOCAL VARIABLES
      DOUBLE COMPLEX DFT(iNODT,iNODT),DFTINV(iNODT,iNODT)
!     DFT = DFT MATRIX
!     DFTINV = INVERSE DFT MATRIX
      DOUBLE COMPLEX SEQ(iNODT),COEF(iNODT)
!     SEQ = LOAD INPUT      
!     COEFF = FOURIER COEFFICIENTS
      DIMENSION FREQ(iNODT)
!     FREQ = FREQUENCY GRID
      DOUBLE COMPLEX TEMP1(iNODT)
!     ===================================================================
      WRITE(*,*) '   ANALYSIS TYPE: DFT'

!     DFT MATRIX
      CALL CDFTMAT(DFT,iNODT)
      CALL CIDFTMAT(DFTINV,iNODT)
      
!     FREQUENCY GRID
      TWID=rDT*iNODT
      CALL FREQGRID(FREQ,TWID,iNODT) 

!     LOAD INPUT
      DO i=1,iNODT
        SEQ(i)=DCMPLX(rLV(i),0.0D0)
      END DO
      
!     FOURIER COEFFCIENTS
      COEF=MATMUL(DFT,SEQ)

!     CALCULATE DISPLACEMENT
      PI=3.141592653589793238462643383279
      PI2=PI*2.0D0
      DO i=1,iNODT
        A=-rMASS*(FREQ(i)*PI2)**2+rSTIFF
        B=(FREQ(i)*PI2)*rDAMPCOEFF
        COEF(i)=(DCMPLX(1.0D0,0.0D0)/DCMPLX(A,B))*COEF(i)
      END DO
      TEMP1=MATMUL(DFTINV,COEF)  
      
!     SAVE DISPLACEMENT 
      DO i=1,iNODT
        rU(i)=DREAL(TEMP1(i))
      END DO     

!     CALCULATE VELOCITY
      DO i=1,iNODT
        COEF(i)=DCMPLX(0.0D0,1.0D0)*COEF(i)*(FREQ(i)*PI2)
      END DO
      TEMP1=MATMUL(DFTINV,COEF)

!     SAVE VELOCITY 
      DO i=1,iNODT
        rUD(i)=DREAL(TEMP1(i))
      END DO      

!     CALCULATE ACCELERATION AND SAVE     
      DO i=1,iNODT
        rUDD(i)=(rLV(i)-rDAMPCOEFF*rUD(i)-rSTIFF*rU(i))/rMASS
      END DO
      
      RETURN
      END SUBROUTINE
!     ===================================================================

!     ===================================================================
      SUBROUTINE ANALYSIS_IMPROVEDDFT(iMAXNOS,rSTIFF,rDT,iNODT,
     & rLV,rWN,rWD,rU,rUD,rUDD,rMASS,rDAMPCOEFF,rDAMPRATIO,iANATYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION
!     * LOCAL VARIABLES
      DOUBLE COMPLEX DFT(iNODT,iNODT),DFTINV(iNODT,iNODT)
!     DFT = DFT MATRIX
!     DFTINV = INVERSE DFT MATRIX
      DOUBLE COMPLEX SEQ(iNODT),COEF(iNODT)
!     SEQ = LOAD INPUT      
!     COEFF = FOURIER COEFFICIENTS
      DIMENSION FREQ(iNODT)
!     FREQ = FREQUENCY GRID
      DOUBLE COMPLEX TEMP1(iNODT)
!     ===================================================================
      WRITE(*,*) '   ANALYSIS TYPE: IMPROVED DFT'

!     SAVE INITIAL CONDITION
      U=rU(1)
      UD=rUD(1)

!     DFT MATRIX
      CALL CDFTMAT(DFT,iNODT)
      CALL CIDFTMAT(DFTINV,iNODT)
      
!     FREQUENCY GRID
      TWID=rDT*iNODT
      CALL FREQGRID(FREQ,TWID,iNODT) 
      
!     LOAD INPUT
      DO i=1,iNODT
        SEQ(i)=DCMPLX(rLV(i),0.0D0)
      END DO
      
!     FOURIER COEFFCIENTS
      COEF=MATMUL(DFT,SEQ)

!     CALCULATE DISPLACEMENT
      PI=3.141592653589793238462643383279
      PI2=PI*2.0D0
      DO i=1,iNODT
        A=-rMASS*(FREQ(i)*PI2)**2+rSTIFF
        B=(FREQ(i)*PI2)*rDAMPCOEFF
        COEF(i)=(DCMPLX(1.0D0,0.0D0)/DCMPLX(A,B))*COEF(i)
      END DO
      TEMP1=MATMUL(DFTINV,COEF)  
      
!     SAVE DISPLACEMENT 
      DO i=1,iNODT
        rU(i)=DREAL(TEMP1(i))
      END DO     

!     CALCULATE VELOCITY
      DO i=1,iNODT
        COEF(i)=DCMPLX(0.0D0,1.0D0)*COEF(i)*(FREQ(i)*PI2)
      END DO
      TEMP1=MATMUL(DFTINV,COEF)

!     CALCULATE HOMOGENEOUS SOLUTION A+B (ADJUSTMENT TERM)
      U=U-rU(1)
      UD=UD-DREAL(TEMP1(1))
      A=U
      B=(UD+rDAMPRATIO*A)/rWD
      
!     ADJUST DISPLACEMENT 
      TIME=0.0D0
      DO i=1,iNODT
        H=DEXP(-rDAMPRATIO*rWN*TIME)*(A*DCOS(rWD*TIME)+B*SIN(rWD*TIME))
        rU(i)=rU(i)+H
        TIME=TIME+rDT
      END DO      

!     SAVE VELOCITY
      TIME=0.0D0
      C=-rDAMPRATIO*rWN*A+B*rWD
      D=-rDAMPRATIO*rWN*B-A*rWD
      DO i=1,iNODT
        H=DEXP(-rDAMPRATIO*rWN*TIME)*(C*DCOS(rWD*TIME)+D*SIN(rWD*TIME))
        rUD(i)=DREAL(TEMP1(i))+H
        TIME=TIME+rDT
      END DO      

!     CALCULATE ACCELERATION AND SAVE     
      DO i=1,iNODT
        rUDD(i)=(rLV(i)-rDAMPCOEFF*rUD(i)-rSTIFF*rU(i))/rMASS
      END DO
      
      RETURN
      END SUBROUTINE
!     ===================================================================

!     ===================================================================
      SUBROUTINE ANALYSIS_STATEFORM(iMAXNOS,rSTIFF,rDT,iNODT,
     & rLV,rWN,rWD,rU,rUD,rUDD,rMASS,rDAMPCOEFF,rDAMPRATIO,iANATYPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION
!     * LOCAL VARIABLES
      DIMENSION A(2,2),B(2,2),ULOAD(2,iNODT),TRANMAT(2,2),TRANMAT2(2,2),
     & XZERO(2)
!     A, B = COEFFICIENT MATRIX OF STATE FORM XD=AX+BU
!     ULOAD = LOAD MATRIX  
!     TRANMAT = STATE-TRANSITION MATRIX
!     XZERO = INITIAL VALUE
      DIMENSION TXH(2),TXP(2)
!     TXH = HOMOGENEOUS SOLUTION OF TIME T
!     TXP = PARTICULAR SOLUTION OF TIME T    
!     ===================================================================
      WRITE(*,*) '   ANALYSIS TYPE: STATEFORM'

!     INITIAL VALUE     
      XZERO(1)=rU(1)
      XZERO(2)=rUD(1)

!     LOAD MATRIX
      ULOAD(:,:)=0.0D0
      DO i=1,iNODT
        ULOAD(1,i)=rLV(i)
      END DO        

!     COEFFICIENT A, B
      A(:,:)=0.0D0
      B(:,:)=0.0D0
      
      A(1,2)=1.0D0
      A(2,1)=-rSTIFF/rMASS      
      A(2,2)=-rDAMPCOEFF/rMASS
      
      B(1,2)=1.0D0
      B(2,1)=1.0D0/rMASS
      
!     INITIAL STATE
      TXH(1)=XZERO(1)
      TXH(2)=XZERO(2)
      TXP(:)=0.0D0
      
!     CALCULATIONS FOR EACH TIME STEP
      DO i=2,iNODT
        CALL STATEFORM_STATE_TRANSITION_MATRIX(TRANMAT,A,rDT)
        TXH=MATMUL(TRANMAT,TXH)
        
        CALL STATEFORM_PARTICULAR(TXP,rDT,A,B,ULOAD,i,iNODT)
     
        rU(i)=TXH(1)+TXP(1)
        rUD(i)=TXH(2)+TXP(2)
     
      END DO

!     CALCULATE ACCELERATION AND SAVE     
      DO i=1,iNODT
        rUDD(i)=(rLV(i)-rDAMPCOEFF*rUD(i)-rSTIFF*rU(i))/rMASS
      END DO      
      
      
      RETURN
      END SUBROUTINE
!     ===================================================================

!     ===================================================================
      SUBROUTINE STATEFORM_STATE_TRANSITION_MATRIX(TRANMAT,A,TIME)
      INCLUDE 'LINK_FNL_SHARED.H'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     * READ
      DIMENSION A(2,2)
!     A = COEFFICIENT MATRIX OF STATE FORM XD=AX+BU
      DOUBLE COMPLEX EVAL(2),EVEC(2,2),TEMP1(2,2),TEMP2(2,2),TEMP3(2,2)
!     EVAL = Complex array of size N containing the eigenvalues of A in decreasing order of magnitude      
!     EVEC = Complex array containing the matrix of eigenvectors. (Output)
!               The J-th eigenvector, corresponding to EVAL(J), is stored in the J-th column. Each
!               vector is normalized to have Euclidean length equal to the value one
!    * RETURN
      DIMENSION TRANMAT(2,2)
!     TRANMAT = STATE-TRANSITION MATRIX
!     ===================================================================
      
!     EIGENVALUES AND EIGENVECTOR
      CALL DEVCRG(2,A,2,EVAL,EVEC,2)      
      
!     STATE-TRANSITION MATRIX
      TEMP1(:,:)=0.0D0
      TEMP2(:,:)=0.0D0
      
      TEMP1(1,1)=CDEXP(EVAL(1)*DCMPLX(TIME,0))
      TEMP1(2,2)=CDEXP(EVAL(2)*DCMPLX(TIME,0))
      
      CALL DLINCG(2,EVEC,2,TEMP2,2)
      TEMP3=MATMUL(MATMUL(EVEC,TEMP1),TEMP2)
      
      TRANMAT(1,1)=DREAL(TEMP3(1,1))
      TRANMAT(1,2)=DREAL(TEMP3(1,2))
      TRANMAT(2,1)=DREAL(TEMP3(2,1))
      TRANMAT(2,2)=DREAL(TEMP3(2,2))
      
      RETURN
      END SUBROUTINE
!     ===================================================================

!     ===================================================================
      SUBROUTINE STATEFORM_PARTICULAR(TXP,rDT,A,B,ULOAD,iNDEX,iNODT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     * READ
      DIMENSION A(2,2),B(2,2),ULOAD(2,iNODT)
!     A, B = COEFFICIENT MATRIX OF STATE FORM XD=AX+BU
!     ULOAD = LOAD MATRIX  
!     * RETURN
      DIMENSION TXP(2)
!     TXP = PARTICULAR SOLUTION OF TIME T  
!     * LOCAL VARIABLE
      DIMENSION TRANMAT(2,2),TEMP1(2,2),rLOAD1(2),rLOAD2(2),
     & TMP(2)
!     TRANMAT = STATE-TRANSITION MATRIX
!     ===================================================================
      
      DO i=1,2
        rLOAD1(i)=ULOAD(i,iNDEX-1)
        rLOAD2(i)=ULOAD(i,iNDEX)
      END DO
      
      CALL STATEFORM_STATE_TRANSITION_MATRIX(TRANMAT,A,-rDT)
      
      TEMP1=MATMUL(TRANMAT,B)
      CALL MATMULTI(TEMP1,2,rLOAD2,2,2)
      CALL MATMULTI(B,2,rLOAD1,2,2)
      
      DO i=1,2
        TMP(i)=(rLOAD2(i)+rLOAD1(i))*rDT/2.0
      END DO
      
      DO i=1,2
        TXP(i)=TXP(i)+TMP(i)
      END DO
      
      CALL STATEFORM_STATE_TRANSITION_MATRIX(TRANMAT,A,rDT)
      
      CALL MATMULTI(TRANMAT,2,TXP,2,2)

      RETURN
      END SUBROUTINE
!     ===================================================================
