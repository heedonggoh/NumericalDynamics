!     =======================================================================
      SUBROUTINE TIMEGRID(TIME,WIDTH,iDISC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TIME(iDISC)
!     =======================================================================
      DX=WIDTH/(iDISC*1.0D0)
      DO i=1,iDISC
        TIME(i)=-WIDTH/2.0D0+WIDTH/(iDISC*1.0D0)+(i-1)*DX
      END DO
                  
      RETURN
      END SUBROUTINE
!     =======================================================================

!     =======================================================================
      SUBROUTINE FREQGRID(FREQ,WIDTH,iDISC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION FREQ(iDISC)
!     =======================================================================
      DX=WIDTH/(iDISC*1.0D0)
      DFREQ=1.0D0/WIDTH
      DO i=1,iDISC
        FREQ(i)=-1.0D0/(2.0D0*DX)+1.0D0/WIDTH+(i-1)*DFREQ
      END DO
                  
      RETURN
      END SUBROUTINE
!     =======================================================================

!     =======================================================================
      SUBROUTINE CDFTMAT(DFTMAT,iORDER)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE COMPLEX DFTMAT(iORDER,iORDER)
      PI=3.141592653589793      
!     =======================================================================
      IF(MOD(iORDER,2).NE.0) THEN
        WRITE(*,*) 'ERROR: SUBROUTINE DFTMAT-ORDER OF THE MATRIXIS ODD!'
        PAUSE
        STOP
      END IF
      
      rORDER=iORDER*1.0D0

      K=-iORDER/2+1
      DO i=1,iORDER
        N=-iORDER/2+1
        DO j=1,iORDER
            DFTMAT(j,i)=DCMPLX(DCOS((2.0D0*PI*N*K)/rORDER)/rORDER,
     &                          -DSIN((2.0D0*PI*N*K)/rORDER)/rORDER)
            N=N+1
        END DO
      K=K+1  
      END DO                    
      RETURN
      END SUBROUTINE
!     =======================================================================

!     =======================================================================
      SUBROUTINE CIDFTMAT(DFTINV,iORDER)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE COMPLEX DFTINV(iORDER,iORDER)
      PI=3.141592653589793      
!     =======================================================================
      IF(MOD(iORDER,2).NE.0) THEN
        WRITE(*,*) 'ERROR: SUBROUTINE DFTMAT-ORDER OF THE MATRIXIS ODD!'
        PAUSE
        STOP
      END IF

      rORDER=iORDER*1.0D0

      K=-iORDER/2+1
      DO i=1,iORDER
        N=-iORDER/2+1
        DO j=1,iORDER
            DFTINV(j,i)=DCMPLX(DCOS((2.0D0*PI*N*K)/rORDER),
     &                          DSIN((2.0D0*PI*N*K)/rORDER))
            N=N+1
        END DO
      K=K+1  
      END DO                    
            
      RETURN
      END SUBROUTINE
!     =======================================================================

