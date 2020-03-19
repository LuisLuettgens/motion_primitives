      SUBROUTINE RBRMDE(IPHASE,TIME,YVEC,NYVEC,PARM,NPARM,
     $    FRHS,NFRHS,IFERR)
C
C         Computes the right hand sides of the robot arm problem
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER    IPHASE,NYVEC,NPARM,NFRHS,IFERR
      DIMENSION  YVEC(NYVEC),PARM(NPARM),FRHS(NFRHS)
c
      parameter (zero=0.d0,one=1.d0,two=2.d0)
c
      parameter (capL = 5.d0)
c
      dimension ucntrl(3)
C
      IFERR = 0
c
c         unload dynamic variables into local variables
c
c         ---algebraic variables
c
      ucntrl(1:3) = yvec(7:9)
c
      phiI = ((capL - yvec(1))**3 + yvec(1)**3)/3.d0
c
      thetaI = phiI*(sin(yvec(5))**2)
c
      frhs(1) = yvec(2)    
c
      frhs(2) = ucntrl(1)/capL
c
      frhs(3) = yvec(4)    
c
      frhs(4) = ucntrl(2)/thetaI
c
      frhs(5) = yvec(6)    
c
      frhs(6) = ucntrl(3)/phiI
C
      RETURN
      END

