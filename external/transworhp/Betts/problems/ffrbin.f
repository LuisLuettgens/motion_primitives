
      SUBROUTINE FFRBIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
c         Free Flying Robot Model of Yoshiyuki Sakawa, Optimal Control 
c         Applications and Methods, 20, 235-248 (1999), suggested by
c         Helmut Maurer.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,ONEEP1=1.D+01,ONEEP3=1.D+03)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(5),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      character  title(3)*60,plbl(mxparm)*80,stlbl(0:mxstat)*80,
     &           clbl(0:mxpcon)*80
c
      INCLUDE '../commons/odeprb.cmn'
c
C     ******************************************************************
C
c
C             Obtain number of grid points for problem.
C
      NPHS = 1
      ICL = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),-2;(la3),-20')
cref      call STMSEQ('(la2),-2;(la3),-3;(la4),-20')
      call STMSEQ('(trp),-2;(la3),-3;(la4),-20')
      NTERM = 0
      NKON = 0
      BIGBND = ONE/HDMCON(5)
c
      TITLE(1) = 'Free Flying Robot, Yoshiyuki Sakawa'
      TITLE(2) = 'Absolute Value Elimination by Slacks'
c
      NPV = 0
      NAV = 4
      NDE = 6
      NCF(1) = NDE
c
C             Initial time and final time fixed.
C
      Y0(0) = zero
      Y1(0) = 12.d0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables (X).
C
      call dfill(nde+nav,zero,y0(1),1)
c
      y0(1) = -10.d0
      y0(2) = -10.d0
      y0(3) = hdmcon(12)/2.d0
      y0(4) = 0.d0
      y0(5) = 0.d0
      y0(6) = 0.d0
C
C             Define terminal conditions for state variables (X).
C
      call dfill(nde+nav,zero,y1(1),1)
c
c     ----default bounds
c
      do 120 j = -1,1
        ylb(j,7) = 0.d0
        ylb(j,8) = 0.d0
        ylb(j,9) = 0.d0
        ylb(j,10) = 0.d0
c
        yub(j,7) = 1.d0
        yub(j,8) = 1.d0
        yub(j,9) = 1.d0
        yub(j,10) = 1.d0
 120  continue
c
c       fix the initial and final states
c
      do 130 j = 1,6
        ylb(-1,j) = y0(j)
        yub(-1,j) = y0(j)
        ylb(1,j) = y1(j)
        yub(1,j) = y1(j)
 130  continue
c
c       ---set up quadrature objective
c
      ncf(3) = ncf(3) + 1
      NTERM = NTERM + 1
      ITERM(1,NTERM) = 0
      ITERM(2,NTERM) = IPHASE
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = -ncf(2) - ncf(3)
      COEF(NTERM) = 1.d0
c
      MAXMIN = -1
C
      RETURN
      END
