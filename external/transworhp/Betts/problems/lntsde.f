

      SUBROUTINE LNTSDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the linear tangent
C         steering system of differential equations.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0)
C
C  Arguments:
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
C
C  Local:
      COMMON /LNTSCM/ ATHRUS,CPI2,ICL
C
C     ******************************************************************
C
C             Compute optimal control.
C
      IF (ICL.EQ.1) THEN
C
C             Compute control and adjoint equations.
C
          DENOM = SQRT(Y(7)**2 + Y(8)**2)
          IF (DENOM.EQ.ZERO)
     &      THEN
              IFERR = 1
              RETURN
            ELSE 
              COSU = -Y(7)/DENOM
              SINU = -Y(8)/DENOM
          ENDIF
          F(5) = ZERO
          F(6) = ZERO
          F(7) = - Y(5)
          F(8) = - Y(6)
      elseif(icl.eq.4) then
c
        beta = atan(p(1) - p(2)*t)
        cosu = cos(beta)
        sinu = sin(beta)
c
      ELSE
          COSU = COS(Y(5))
          SINU = SIN(Y(5))
      ENDIF
C
C             Compute state equations.
C
      F(1) = Y(3)
      F(2) = Y(4)
      F(3) = ATHRUS*COSU
      F(4) = ATHRUS*SINU
c
c       nonzero gravity approximation:
c       F(4) = ATHRUS*SINU - 32.2d0*(3500.d0/(Y(2)+3500.d0))**2
c
      if(icl.eq.3) then
        if(t.gt..2d0) f(2) = -10.d0*f(2) 
      endif
C
      IFERR = 0
C
      RETURN
      END
