!     ===================================================================
!     NUMERICAL EVALUATION OF DYNAMIC RESPONSE
!     CREATED 2010.11.18 UPDATED 2011.05.16
!     HEEDONG GOH <wellposed@gmail.com>
!     * ANALYSIS TYPE *
!      1. LINEAR INTERPOLATION
!      2. CENTRAL DIFFERENCE METHOD
!      3. AVERAGE ACCELERATION NEWMARK'S METHOD
!      4. LINEAR ACCELERATION NEWMARK'S METHOD
!      5. DFT METHOD
!      6. IMPROVED DFT METHOD
!      7. STATE FORM METHOD
!     * LOAD TYPE *
!      1. HARMONIC LOADING
!      2. PULSE LOADING
!      3. EARTHQUAKE LOADING(EL CENTRO 1940 NORTH SOUTH COMPONENT)
!     MAIN PROGRAMME
      PROGRAM NUMDYNAMICS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     ===================================================================
      PARAMETER(iMAXNOS=20000,rPI=3.14159265358979323846264338327950)
!     iMAXNOS = MAX NUMBER OF TIME STEPS
!     rPI = pi
      DIMENSION rLV(iMAXNOS),rU(iMAXNOS),rUD(iMAXNOS),rUDD(iMAXNOS),
     & rUDDG(iMAXNOS)
!     rLV = LOAD VECTOR 
!     rU = DISPLACEMENT
!     rUD = VELOCITY
!     rUDD = ACCELERATION    
!     rUDDG = GROUND MOTION
!     ===================================================================
      
      WRITE(*,*)
      WRITE(*,*) '# PROGRAMME: NUMERICAL EVALUATION OF DYNAMIC RESPONSE'
      WRITE(*,*)
      
      CALL INPUT(iMAXNOS,iANATYPE,rMASS,rDAMPCOEFF,rSTIFF,rDT,
     & iNODT,iLOADTYPE,rPZERO,rFREQ,rU,rUD,iPULSETYPE,rDURA,rUDDG)
       
      CALL PREPROCESS(iMAXNOS,rPI,rMASS,rDAMPCOEFF,rSTIFF,rDT,
     & iNODT,iLOADTYPE,rPZERO,rFREQ,rUDDG,
     & rLV,rWN,rDAMPRATIO,rWD,rTN,rTD,rDURA,iPULSETYPE)
     
      CALL ANALYSIS(iMAXNOS,rSTIFF,rDT,iNODT,iLOADTYPE,iANATYPE,
     & rLV,rWN,rDAMPRATIO,rWD,rU,rUD,rUDD,rMASS,rDAMPCOEFF)
     
      CALL OUTPUT(iMAXNOS,iANATYPE,iLOADTYPE,iNODT,rDT,rPI,
     & rMASS,rDAMPCOEFF,rSTIFF,rPZERO,rFREQ,rDAMPRATIO,rWN,rWD,rTN,rTD,
     & rLV,rU,rUD,rUDD,rUDDD,rDURA,iPULSETYPE)
     
      WRITE(*,*)
      WRITE(*,*) '# ANALYSIS COMPLETED'
      WRITE(*,*) 
      
      END PROGRAM

!     ===================================================================
      
