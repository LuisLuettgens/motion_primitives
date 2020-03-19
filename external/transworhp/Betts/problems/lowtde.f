

      SUBROUTINE LOWTDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the low thrust orbit
c         transfer differential equations.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,one=1.d0)
C
C  Arguments:
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
C
C  Local:
      COMMON /LOWTCM/ ATHRUS,CPI2,ICL
C
C     ******************************************************************
C
C             Compute optimal control.
C
      radius = y(1)
      theta = y(2)
      vsubr = y(3)
      vsubt = y(4)
      COSBET = COS(Y(5))
      SINBet = SIN(Y(5))
C
C             Compute state equations.
C
      rdot = vsubr
      thetdt = vsubt/radius
      vrdot = vsubt**2/radius - one/(radius**2) + athrus*sinbet
      vtdot = - vsubt*vsubr/radius + athrus*cosbet
c
      F(1) = rdot
      F(2) = thetdt
      F(3) = vrdot
      F(4) = vtdot
C
      IFERR = 0
C
      RETURN
      END
