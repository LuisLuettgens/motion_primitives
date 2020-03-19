      SUBROUTINE RBRMIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the robot arm problem  ------  
c             Ref: Thesis of Monika Mossner-Beigel (Heidelberg Univ).
c             Example 8 (p. 18) in "Benchmarking Optimization Software 
c             with COPS 3.0".  Mathematics and Computer Science Division,  
c             Argonne National Laboratory,  Technical Report ANL/MCS-TM-273, 
c             February 2004.              
c
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,two=2.D0)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT(2),MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           clbl(0:mxpcon)*80
c
      INCLUDE '../commons/odeprb.cmn'
C
C     ******************************************************************
c
      iclslb = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
c
      picons = hdmcon(12)
c
      TITLE(1) = 'Robot Arm Problem'
      TITLE(2) = 'Minimum Time Motion'
c
      NDE = 0
      NPV = 0
      NAV = 0
C
      NTERM = 0
      NKON = 0
      BIGBND = ONE/HDMCON(5)
c
C     ==================================================================
C     ========Differential Variables====================================
C     ==================================================================
c
      STLBL(nde) = 'TIME    Time'
c
      y0(nde) = zero
      y1(nde) = one  
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'Y1      State y1'
c
      y0(nde) = 9.d0/two
      y1(nde) = 9.d0/two
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'Y2      State y2'
c
      y0(nde) = zero 
      y1(nde) = zero 
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'Y3      State y3'
c
      y0(nde) = zero 
      y1(nde) = (two*picons)/3.d0
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'Y4      State y4'
c
      y0(nde) = zero 
      y1(nde) = zero 
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'Y5      State y5'
c
      y0(nde) = picons/4.d0
      y1(nde) = picons/4.d0
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'Y6      State y6'
c
      y0(nde) = zero 
      y1(nde) = zero  
c
c     ------------------------------
c
      NCF(1) = NDE
c
c     ==================================================================
c     ========Algebraic Variables=======================================
c     ==================================================================
c
      nav = nav + 1
      STLBL(nde+nav) = 'U1      Control 1'
c
      y0(nde+nav) = zero
      y1(nde+nav) = zero
c
      ylb(-1:1,nde+nav) = -one
      yub(-1:1,nde+nav) = +one
c
c     ------------------------------
c
      nav = nav + 1
      STLBL(nde+nav) = 'U2      Control 2'
c
      y0(nde+nav) = zero
      y1(nde+nav) = zero
c
      ylb(-1:1,nde+nav) = -one
      yub(-1:1,nde+nav) = +one
c
c     ------------------------------
c
      nav = nav + 1
      STLBL(nde+nav) = 'U3      Control 3'
c
      y0(nde+nav) = zero
      y1(nde+nav) = zero
c
      ylb(-1:1,nde+nav) = -one
      yub(-1:1,nde+nav) = +one
c
c     ==================================================================
C     ========Boundary Conditions=======================================
c     ==================================================================
c
c          fixed initial time and state
c
      ylb(-1,0:nde) = y0(0:nde)
      yub(-1,0:nde) = y0(0:nde)
c
c          fixed final state
c
      ylb(+1,1:nde) = y1(1:nde)
      yub(+1,1:nde) = y1(1:nde)
c
c     ==================================================================
C     ========Objective Function========================================
C     ==================================================================
C
c
      CLBL(0) = 'TFINAL  Final Time'
C
      maxmin = -1
c
      NTERM = NTERM + 1
      ITERM(1,NTERM) = 0
      ITERM(2,NTERM) = 1
      ITERM(3,NTERM) = 1
      ITERM(4,NTERM) = 0    
      COEF(NTERM) = one
C
C ========================================================================
c
c
  
      RETURN
      END
