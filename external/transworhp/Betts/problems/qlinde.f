

      SUBROUTINE QLINDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the quadratic linear
C         controller system of differential equations.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,HALF=0.5D0)
C
C  Arguments:
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(*),F(NF)
C
      COMMON /QLCOM/ FMTRX(6,6),GMTRX(6,3),AMTRX(6,6),BMTRX(3,3),
     &               SFMTRX(6,6),ICL,NDEQ,NALG
C
C     ******************************************************************
C
      IFERR = 0
C
C             Compute F = FMTRX*Y + GMTRX*U.
C
      IF (ICL.EQ.1) THEN
c
C             Quadrature = 0.5*( (Y**t)*AMTRX*Y + (U**t)*BMTRX*U) )
C
          CALL HDMVPS(1,NDEQ,NDEQ,AMTRX,6,Y,F,IFERR)
          QF = DDOT(NDEQ,Y,1,F,1)
          CALL HDMVPS(1,NALG,NALG,BMTRX,3,Y(NDEQ+1),F,IFERR)
          QF = HALF*(QF+DDOT(NALG,Y(NDEQ+1),1,F,1))
C
          CALL HDMVPS(1,NDEQ,NDEQ,FMTRX,6,Y,F,IFERR)
          CALL HDMVPS(2,NDEQ,NALG,GMTRX,6,Y(NDEQ+1),F,IFERR)
          F(NDEQ+1) = QF
c
      ELSEIF(ICL.EQ.2) THEN
c
          NS = NDEQ - 1
C
C             Control integral:  F(7) = 0.5*(U**t)*BMTRX*U
C
          CALL HDMVPS(1,NALG,NALG,BMTRX,3,Y(NDEQ+1),F,IFERR)
          F(NDEQ) = HALF*DDOT(NALG,Y(NDEQ+1),1,F,1)
C
          CALL HDMVPS(1,NS,NS,FMTRX,6,Y,F,IFERR)
          CALL HDMVPS(2,NS,NALG,GMTRX,6,Y(NDEQ+1),F,IFERR)
C
      ELSEIF (ICL.EQ.3) THEN
C
C             Quadrature = 0.5*( (Y**t)*AMTRX*Y + (U**t)*BMTRX*U) )
C
          CALL HDMVPS(1,NDEQ,NDEQ,AMTRX,6,Y,F,IFERR)
          QF = DDOT(NDEQ,Y,1,F,1)
          CALL HDMVPS(1,NALG,NALG,BMTRX,3,Y(NDEQ+1),F,IFERR)
          QF = HALF*(QF+DDOT(NALG,Y(NDEQ+1),1,F,1))
C
          CALL HDMVPS(1,NDEQ,NDEQ,FMTRX,6,Y,F,IFERR)
          CALL HDMVPS(2,NDEQ,NALG,GMTRX,6,Y(NDEQ+1),F,IFERR)
          F(NDEQ+1) = .1D0*Y(1) + .2D0*Y(2)
          F(NDEQ+2) = QF
C
      ELSEIF (ICL.EQ.4) THEN
C
C             QUADRATURE OBJECTIVE -- NO ODE'S OR ALG. EQS.
C
          PI = HDMCON(12)
          UDES = SIN(2.D0*PI*T)
          F(1) = (UDES - Y(1))**2
C
      ENDIF
C
      RETURN
      END
