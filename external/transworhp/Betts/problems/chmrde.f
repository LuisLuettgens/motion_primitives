

      SUBROUTINE CHMRDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the chemical reator
c         differential equations.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,one=1.d0)
C
C  Arguments:
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
C
C  Local:
      COMMON /CHMRCM/ alwr,aupr,tfinal,xzero,yzero,rhocof,expk
C
C     ******************************************************************
c
      iferr = 0
C
C             Compute optimal control.
C
      xprod = y(1)
      yprod = y(2)
      arate = y(3)
c
      if(arate.le.zero) then
        iferr = 1
        return
      endif
c
      bcoef = rhocof*(arate**expk)
C
C             Compute state equations.
C
      xdot = -arate*xprod
      ydot = arate*xprod - bcoef*yprod
c
      F(1) = xdot
      F(2) = ydot
C
      IFERR = 0
C
      RETURN
      END
