
      SUBROUTINE TRANIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C             TRAIN EXAMPLE PROBLEM
C
C             REF:  Vanderbei, "Case Studies in Trajectory Optimization:
c                   Trains, Planes, and other Pastimes", Operations Research
c                   and Financial Engineering, Princeton Univ., July 27, 2000
C
C     ******************************************************************
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Arguments:
C
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           CLBL(0:MXPCON)*80
c
      INCLUDE '../commons/odeprb.cmn'
C
C     ******************************************************************
c
      bigbnd = 1.d0/hdmcon(5)
      NTERM = 0
      NKON = 0
C
C     ----PROBLEM TITLES AND LABELS (OPTIONAL)
C
      TITLE(1) = 'Train Problem'
      TITLE(2) = 'Minimum Fuel Cost--With Regularization'
      STLBL(0) = 'TIME    Time'
      STLBL(1) = 'XPOS    Position Along the Track'
      STLBL(2) = 'VPOS    Velocity Along the Track'
      STLBL(3) = 'USUBA   Acceleration from Engines'
      STLBL(4) = 'USUBB   Deceleration from Brakes'
c
C     ----NUMBER OF DIFFERENTIAL EQUATIONS (STATES)   
C
      NDE = 2
C
C     ----NUMBER OF ALGEBRAIC VARIABLES (CONTROLS)   
C
      NAV = 2
      NCF(1) = NDE
C
C     ----GUESS FOR INITIAL TIME AND BOUNDARY CONDITION
C
      Y0(0) = 0.D0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
C
C     ----GUESS FOR FINAL TIME AND BOUNDARY CONDITION
C
      T1 = 4.8D0
      Y1(0) = T1
      YLB(1,0) = t1
      YuB(1,0) = t1
C
C     ----DEFINE INITIAL CONDITIONS FOR STATE VARIABLES, Y.
C
      Y0(1) = 0.D0
      Y0(2) = 0.D0
C
C     ----FIX THE INITIAL STATES
C
      DO 110 I=1,nde
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
  110 CONTINUE
C
C     ----DEFINE FINAL CONDITIONS FOR STATE VARIABLES, Y.
C
      Y1(1) = 6.D0
      Y1(2) = 0.D0
C
C     ----DEFINE INITIAL AND FINAL CONTROL ANGLES.
C
      Y0(3) = 1.D0
      Y0(4) = 1.D0
      Y1(3) = 1.D0
      Y1(4) = 1.D0
C
C     ----FIX THE FINAL VALUES FOR STATES
C
      DO 120 I=1,nde
        YLB(1,I) = Y1(I) 
        YUB(1,I) = Y1(I)
  120 CONTINUE
c
c     ----Bound the controls
c
      do 130 jj = -1,1
        YLB(jj,3) = 0.d0
        YUB(jj,3) = 10.d0
        YLB(jj,4) = 0.d0 
        YUB(jj,4) = 2.d0
 130  continue
C
C     ----DEFINE work AS THE OBJECTIVE TO BE MINIMIZED
C
      ncf(3) = ncf(3) + 1
      NTERM = NTERM + 1
      ITERM(1,nterm) = 0
      ITERM(2,nterm) = iphase
      ITERM(3,nterm) = 0
      ITERM(4,nterm) = -ncf(2)-ncf(3)
      COEF(nterm) = 1.D0
C
      MAXMIN = -1
C
      RETURN
      END
