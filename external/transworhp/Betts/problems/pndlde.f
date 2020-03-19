

      SUBROUTINE PNDLDE(IPHASE,T,Y,NY,P,NP,Frhs,NFrhs,IFERR)
C
C         COMPUTES THE RIGHT HAND SIDES OF THE pendulum
C         DIFFERENTIAL EQUATIONS.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ARGUMENTS:
C
      INTEGER    IPHASE,NY,NP,NFrhs,IFERR
      DIMENSION  Y(NY),P(NP),Frhs(NFrhs)
      COMMON  /pndlCM/ ICLSS
c
      parameter (one=1.d0,two=2.d0,gcon=9.81d0)
C
C     ******************************************************************
C
C     ----SET FUNCTION ERROR FLAG.
C
      IFERR = 0
c
      u = y(6)
c
      frhs(1) = y(3)
      frhs(2) = y(4)
      frhs(3) = - two*y(5)*y(1) + u*y(2)
      frhs(4) = - gcon - two*y(5)*y(2) - u*y(1)
c
      if(iclss.eq.1) then
c
c         path constraint
c
        frhs(5) = y(3)**2 + y(4)**2 - two*y(5) - gcon*y(2)
c
      else
c
        frhs(5) = y(3)*frhs(3) + y(4)*frhs(4) - gcon*frhs(2)/two
c
      endif
c
c         quadrature function
c
      frhs(6) = u**2
c
      RETURN
      END
