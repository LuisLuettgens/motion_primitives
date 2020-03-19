      
      SUBROUTINE AQUAIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
c         Underwater vehicle problem from "SQP-methods for solving optimal
c         control problems with control and state constraints: adjoint
c         variables, sensitivity analysis and real-time control,"
c         Christof Buskens, Helmut Maurer,  JCAM, pp 85-108,
c         Vol 120, Aug. 2000.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,ONEEP1=1.D+01,ONEEP3=1.D+03)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      character  title(3)*60,plbl(mxparm)*80,stlbl(0:mxstat)*80,
     &           clbl(0:mxpcon)*80
      common /aquacm/ iclass
c
      INCLUDE '../commons/odeprb.cmn'
C
C     ******************************************************************
C
C             Obtain number of grid points for problem.
C
      pi = hdmcon(12)
c
      NPHS = 1
      ngrid = itpgrd(iprob)
      method = itpmet(iprob)
      nstg = itpstg(iprob)
      iclass = itpcls(iprob)
      NTERM = 0
      NKON = 0
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),3;(la4),20')
      call STMSEQ('(trp),2;(la3),3;(la4),20')
c
      TITLE(1) = 'Underwater Vehicle, Buskens, Maurer'
      if(iclass.eq.1) then
        TITLE(2) = 'Minimum Energy'
      elseif(iclass.eq.2) then
        TITLE(2) = 'Minimum Energy--Sum of Quadrature Terms'
      endif
c
      NPV = 0
      NAV = 4
      NDE = 10
      NCF(1) = NDE
c
C             Initial time and final time fixed.
C
      Y0(0) = zero
      Y1(0) = one
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables (X).
C
      call dfill(nde+nav,zero,y0(1),1)
      y0(3) = .2d0
      y0(4) = pi/2.d0
      y0(5) = .1d0
      y0(6) = -pi/4.d0
      y0(7) = 1.d0
      y0(9) = .5d0
      y0(10) = .1d0
C
C             Define terminal conditions for state variables (X).
C
      call dfill(nde+nav,zero,y1(1),1)
      y1(1) = 1.d0
      y1(2) = .5d0
      y1(4) = pi/2.d0
c
c     ----default bounds
c
      do 120 j = -1,1
        ylb(j,4) = -.02d0 + pi/2.d0
        yub(j,4) =.02d0 + pi/2.d0
        do 110 i = nde+1,nde+nav
          ylb(j,i) = -15.d0
          yub(j,i) = 15.d0
 110    continue
 120  continue
c
c       fix the initial and final states
c
      do 130 j = 1,nde
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
      ITERM(4,NTERM) = -ncf(3)
      COEF(NTERM) = 1.d0
c
      if(iclass.ne.1) then
        do 140 ii = 2,4
          ncf(3) = ncf(3) + 1
          NTERM = NTERM + 1
          ITERM(1,NTERM) = 0
          ITERM(2,NTERM) = IPHASE
          ITERM(3,NTERM) = 0
          ITERM(4,NTERM) = -ncf(3)
          COEF(NTERM) = 1.d0
 140    continue
      endif
c
      MAXMIN = -1
C
      RETURN
      END
