


      SUBROUTINE LOWTIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the low-thrust spacecraft orbit transfer 
c             problem
C
C             REF:  This is the example problem in "Direct Optimization
c                   Using Collocation Based on High-Order Gauss-Lobatto
c                   Quadrature Rules," by Albert Herman, and Bruce Conway
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FIVE=5.D0)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           clbl(0:mxpcon)*80
C
C  Local:
      COMMON /LOWTCM/ ATHRUS,CPI2,ICL
C
      INCLUDE '../commons/odeprb.cmn'
C
C     ******************************************************************
C
C             Obtain number of grid points for problem.
C
      NPHS = 1
      ICL = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),2;(la3),3;(la4),20')
cref      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
      call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
c
      TITLE(1) = 'Low Thrust Orbit Transfer, Herman/Conway'
      TITLE(2) = 'Three Rev. Solution'
C
      STLBL(0) = 'TIME    Time'
      STLBL(1) = 'RADIUS  Radius'
      STLBL(2) = 'THETA   Range Angle'
      STLBL(3) = 'VSUBR   Radial Velocity'
      STLBL(4) = 'VSUBT   Tangential Velocity'
      STLBL(5) = 'BETA    Thrust Angle (rad)'
c
      NPV = 0
      INIT = 1
      NDE = 4
      NAV = 1
      NCF(1) = NDE
      npf(1) = 0
      npf(2) = 1
      CR2D = HDMCON(16)
      pi = hdmcon(12)
      CPI2 = HDMCON(12)/TWO
      ATHRUS = .01D0
C
C             Default State and Control Bounds
C
      DO 75 I=-1,1
c
c           ---radius
c
        YLB(I,1) = .5d0
        YUB(I,1) = 5.d0
c
c           ---range angle
c
        YLB(I,2) = zero
        YUB(I,2) =  4.d0*(two*pi)
c
c           ---radial velocity
c
        YLB(I,3) = -10.d0
        YUB(I,3) = 10.d0
c
c           ---tangential velocity
c
        YLB(I,4) = zero
        YUB(I,4) = 10.d0
c
c           ---thrust angle
c
        YLB(I,5) = -CPI2
        YUB(I,5) =  CPI2
c
 75   CONTINUE
C
C             Initial time and final time (fixed).
C
      Y0(0) = ZERO
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
c
      T1 = 50.d0
      Y1(0) = T1
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables, Y.
c
C                 Y(1) = Radius
C                 Y(2) = Range Angle
C                 Y(3) = Radial Velocity
C                 Y(4) = Tangential Velocity
C
C             Declare initial values as fixed.
C
      Y0(1) = 1.1d0
      Y0(2) = ZERO
      Y0(3) = ZERO
      Y0(4) = one/sqrt(1.1d0)
      DO 20 I=1,4
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
   20 CONTINUE
C
C             Define terminal conditions for state variables.
C
      call dcopy(4,y0(1),1,y1(1),1)
C
C             Define initial and final control angles.
C
      Y0(5) = zero
      Y1(5) = zero
C
C             Define final specific energy as objective.
C
      ITERM(1,1) = 0
      ITERM(2,1) = 1
      ITERM(3,1) = 1
      ITERM(4,1) = -1
      COEF(1) = ONE
C
      MAXMIN = -1
C
C             Load state scale weights for output.
C
      DO 150 I=0,NDE
        STSKL(I,1) = ONE
  150 CONTINUE
      DO 152 I=NDE+1,NDE+NAV
        STSKL(I,1) = one
  152 CONTINUE
C
      RETURN
      END
