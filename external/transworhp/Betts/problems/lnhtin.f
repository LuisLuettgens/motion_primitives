
      SUBROUTINE LNHTIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
c         Optimal Control of an HIV Chemotherapy Model (Kirschner/Lenhart/Serbin)
c               D. Kirschner, S. Lenhart, S. Serbin:
c               Optimal control of the chemotherapy of HIV,
c               J. of Mathematical Biology, Vol. 35 (1997), 775-792 .
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
      common /lnhtcm/ icls
c
C     ******************************************************************
C
c
C             Obtain number of grid points for problem.
C
      NPHS = 1
      ICLs = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
c
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),20')
      call STMSEQ('(trp),2;(la3),20')
      NTERM = 0
      NKON = 0
      BIGBND = ONE/HDMCON(5)
c
      TITLE(1) = 'HIV Chemotherapy Model'
      TITLE(2) = 'Kirschner/Lenhart/Serbin'
c
      NPV = 0
      NAV = 1
      if(icls.eq.1) then
        NDE = 5
      else
        nde = 4
      endif
      NCF(1) = NDE
c
C             Initial time and final time fixed.
C
      Y0(0) = zero
      Y1(0) = 500.d0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables (X).
C
      call dfill(nde+nav,zero,y0(1),1)
c
      Y0(1) = 9.82d2
      Y0(2) = 0.05d0
      Y0(3) = 6.2d-4
      Y0(4) = 0.02d0
C
C             Define terminal conditions for state variables (X).
C
      call dcopy(nde+nav,y0(1),1,y1(1),1)
c
c     ----default bounds
c
      do 120 j = -1,1
        ylb(j,nde+1) = 0.d0
        yub(j,nde+1) = 1.d0
c
 120  continue
c
c       fix the initial states
c
      do 130 j = 1,nde
        ylb(-1,j) = y0(j)
        yub(-1,j) = y0(j)
 130  continue
c
      if(icls.eq.1) then
c
c       ---set up objective
c
        NTERM = NTERM + 1
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = IPHASE
        ITERM(3,NTERM) = 1
        ITERM(4,NTERM) = nde
        COEF(NTERM) = 1.d0
c
      elseif(icls.eq.2) then
c
c       ---set up objective
c
        ncf(3) = ncf(3) + 1
        NTERM = NTERM + 1
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = IPHASE
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = - ncf(3)
        COEF(NTERM) = 1.d-5
c
      endif
c
      MAXMIN = -1
C
      RETURN
      END
