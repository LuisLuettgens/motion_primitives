
      SUBROUTINE ASHRIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the Ascher Problem, Example 9.2, p. 371
C
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
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
C
C  Local:
C
      INCLUDE '../commons/odeprb.cmn'
C
      COMMON /ASHRCM/ EPS,ICL,NDEQ,NGPT
C
C     ******************************************************************
C
C             Obtain number of grid points for problem.
C
      NPHS = 1
      icl = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
      eps = tpseed(iprob)
c
c
      if(icl.lt.3) then
        TITLE(1) = 'Ascher''s equation, Ex. 9.2, p. 371'
cold      call STMSEQ('(la2),2;(la3),20')
cold        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref        call STMSEQ('(la2),2;(la3),3;(la4),20')
        call STMSEQ('(trp),2;(la3),3;(la4),20')
      elseif(icl.eq.3) then
        TITLE(1) = 'Ascher''s equation, Ex. 10.4, p. 394'
cold      call STMSEQ('(la2),2;(la3),20')
cold        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),3;(la4),20')
        call STMSEQ('(trp),2;(la3),3;(la4),20')
      elseif(icl.eq.4) then
        TITLE(1) = 'Stiff ODE, BCSLIB (HDSODE Example)'
cold      call STMSEQ('(la2),2;(la3),3;(la4),20')
cref      call STMSEQ('(la2),2;(la3),20')
        call STMSEQ('(trp),2;(la3),20')
      elseif(icl.eq.5) then
        TITLE(1) = 'Brusselator, Hairer, Norsett, Wanner, p. 170'
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),3;(la4),20')
        call STMSEQ('(trp),2;(la3),3;(la4),20')
      elseif(icl.eq.6) then
        TITLE(1) = 'Brusselator, Hairer, Norsett, Wanner, p. 170'
        TITLE(2) = 'Slack Variable Formulation'
cold      call STMSEQ('(la2),2;(la3),20')
cold        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),3;(la4),20')
        call STMSEQ('(trp),2;(la3),3;(la4),20')
      endif
c
      NPV = 0
      NAV = 0
      NDE = 2
      NCF(1) = NDE
      nterm = 0
c
      if(icl.eq.6) then
        nav = 4
      endif
c
c         default bounds
c
      do 120 kk = -1,1
        ylb(kk,1) = -5.d0
        yub(kk,1) =  5.d0
        if(icl.eq.5) then
          ylb(kk,1) = -10.d0
          yub(kk,1) =  10.d0
          ylb(kk,2) = -10.d0
          yub(kk,2) =  10.d0
        elseif(icl.eq.6) then
          ylb(kk,1) = -10.d0
          yub(kk,1) =  10.d0
          ylb(kk,2) = -10.d0
          yub(kk,2) =  10.d0
          do 115 jj = 3,6
            ylb(kk,jj) = 0.d0
            ylb(1,jj) = 0.d0
            yub(1,jj) = 0.d0
 115      continue
        else
          ylb(kk,2) = -2500.d0
          yub(kk,2) =  2500.d0
        endif
 120  continue
C
C             Initial time and final time fixed.
C
      if(icl.le.3) then
        Y0(0) = -one
        Y1(0) = ONE
      elseif(icl.eq.4) then
        Y0(0) = zero
        Y1(0) = 5.d0
      elseif(icl.ge.5) then
        Y0(0) = zero
        Y1(0) = 20.d0
      endif
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables (X).
C
      Y0(1) = -2.d0
      Y0(2) = zero
      if(icl.eq.4) then
        Y0(1) = zero
        Y0(2) = one
      elseif(icl.eq.5) then
        Y0(1) = 1.5d0
        Y0(2) = 3.d0
      elseif(icl.eq.6) then
        Y0(1) = 1.5d0
        Y0(2) = 3.d0
        call dfill(nav,0.d0,y0(3),1)
        call dfill(nav,0.d0,y1(3),1)
      endif
C
C             Define terminal conditions for state variables (X).
C
      Y1(1) = ZERO
      Y1(2) = zero
c
      if(icl.eq.3) then
        y0(1) = -one
        y1(1) = one
      endif
C
C             Fix 1st initial state.
C
      YLB(-1,1) = Y0(1)
      YUB(-1,1) = Y0(1)
c
      if(icl.eq.2) then
c
        TITLE(2) = 'Fixed initial slope y2'
C
C             Fix 2nd initial state.
C
        YLB(-1,2) = y0(2)
        YUB(-1,2) = y0(2)
c
      elseif(icl.eq.1.or.icl.eq.3) then
c
        TITLE(2) = 'Fixed initial and final y1'
C
C             Fix 1st final state.
C
        YLB(1,1) = Y1(1)
        YUB(1,1) = Y1(1)
c
      elseif(icl.ge.4) then
c
        YLB(-1,2) = y0(2)
        YUB(-1,2) = y0(2)
c
      endif
C
      MAXMIN = 0
C
      init = 1
c
      if(icl.eq.6) then
C
C         Slack deviation integral
C
        NCF(3) = NCF(3) + 1
        NTERM = NTERM + 1
        CLBL(0) = 'SLKERR  Slack Variable Deviation integral'
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = -1
c
        maxmin = -1
c
      endif
C
C             Load state scale weights for output.
C
      DO 150 I=1,1+NDE+NAV
        STSKL(I,1) = ONE
  150 CONTINUE
C
      IER = 0
C
      RETURN
      END
