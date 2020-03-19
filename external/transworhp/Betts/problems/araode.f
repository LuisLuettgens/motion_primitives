

      SUBROUTINE ARAODE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         COMPUTES THE RIGHT HAND SIDES OF THE A. RAO
C         SYSTEM OF DIFFERENTIAL EQUATIONS.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ARGUMENTS:
C
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
      COMMON  /ARAOCM/ TFINAL,ICLSS,NGRDGS
C
C     ******************************************************************
C
C     ----SET FUNCTION ERROR FLAG.
C
      IFERR = 0
C
      IF(ICLSS.LE.2) THEN
C
C     ----COMPUTE STATE EQUATIONS.
C
        F(1) = -Y(1)**3 + Y(NY)
C
C     ----COMPUTE QUADRATURE EQUATIONS.
C
        F(2) = Y(1)**2 + Y(NY)**2 
C
      ELSEIF(ICLSS.EQ.3) THEN
C
        UCNT = - Y(2)/2.D0
C
C     ----COMPUTE STATE EQUATIONS.
C
        F(1) = -Y(1)**3 + UCNT
C
C     ----COMPUTE ADJOINT EQUATION
C
        F(2) = -2.D0*Y(1) + 3.D0*Y(2)*Y(1)**2
C
      ENDIF
C
      RETURN
      END
