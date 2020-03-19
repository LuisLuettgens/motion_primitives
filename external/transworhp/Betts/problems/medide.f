

      SUBROUTINE MEDIDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         COMPUTES THE RIGHT HAND SIDES OF THE LINEAR TANGENT
C         STEERING SYSTEM OF DIFFERENTIAL EQUATIONS.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ARGUMENTS:
C
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
c
      common /medicm/ smlL,capJ,ngrtot,ipmedi
C
C     ******************************************************************
C
C     ----SET FUNCTION ERROR FLAG.
C
      IFERR = 0
C
      x = y(1)
      v = y(2)
      u = y(3)
c
      f(1) = v
      f(2) = u
c
      if(ipmedi.lt.0) then
c
c         path constraint
c
        f(3) = x
c
c         quadrature objective
c
        f(4) = u**2
c
      else
c
c         quadrature objective
c
        f(3) = u**2
c
      endif
C
      RETURN
      END
