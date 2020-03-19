


      SUBROUTINE BRGRIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the Burgers' equation example.
C
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,ONEEP1=1.D+01,ONEEP3=1.D+03)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT(2),MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
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
      COMMON /BRGRCM/ EPS,ICL,NDEQ,NGPT
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
      EPS = TPSEED(IPROB)
c
cold      call STMSEQ('(la2),2;(la3),20')
cref      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
      call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      TITLE(1) = 'Burgers'' equation'
c
      if(icl.eq.1) then
        TITLE(2) = 'Boundary layer'
        NPV = 0
        NAV = 0
        NDE = 2
        NCF(1) = NDE
        NDEQ = NDE
        NGPT = NGRID
      elseif(icl.eq.2) then
        title(2) = 'Equidistributed Error Formulation'
        NPV = 1
        NAV = 0
        NDE = 3
        NCF(1) = NDE
        NDEQ = NDE
        NGPT = NGRID
      elseif(icl.eq.3) then
        title(2) = 'Elastic Programming Formulation'
        NPV = 0
        NAV = 4
        NDE = 2
        NCF(1) = NDE
        ncf(3) = 1
        NDEQ = NDE
        NGPT = NGRID
      elseif(icl.eq.4) then
        TITLE(2) = 'Optimal Grid Distribution'
        NPV = 0
        NAV = 0
        NDE = 2
        NCF(1) = NDE
        NDEQ = NDE
        NGPT = 10
        ngrid = ngpt
        method = 3
        nstg = 1
      endif
C
C             Initial time and final time fixed.
C
      Y0(0) = ZERO
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      Y1(0) = ONE
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables (X).
C
      REPS = ONE/EPS
      R2EPS = 2.D0/EPS
c      EXPR = EXP(REPS)
c      EXPMR = ONE/EXPR
c      HYPTAN = (EXPR-EXPMR)/(EXPR+EXPMR)
      expr2 = exp(max(hdmcon(10),-r2eps))
      hyptan = (one - expr2)/(one + expr2)
      Y0(1) = 2.D0*HYPTAN
      Y0(2) = -R2EPS*(ONE-HYPTAN**2)
C
C             Fix 1st initial state.
C
      YLB(-1,1) = Y0(1)
      YUB(-1,1) = Y0(1)
C
C             Set lower bound of 0 on first state.
C             (The constraint solve is almost impossible without this)
C
      YLB(0,1) = ZERO
C
C             Define terminal conditions for state variables (X).
C
      Y1(1) = ZERO
      Y1(2) = -R2EPS
C
C             Fix 1st final state.
C
      YLB(1,1) = Y1(1)
      YUB(1,1) = Y1(1)
C
      MAXMIN = 0
C
      if(icl.eq.1) then
        init(1) = 1
      elseif(icl.eq.2) then
c
c           fix the last state variable (time) at initial and final pts.
c
        y0(3) = y0(0)
        y1(3) = y1(0)
        YLB(-1,3) = Y0(3)
        YUB(-1,3) = Y0(3)
        YLB(1,3) = Y1(3)
        YUB(1,3) = Y1(3)
c
c           initialize the parameter
c
        p0(1) = one/(y1(0)-y0(0))
c
        init(1) = 1
c
      elseif(icl.eq.3) then
c
c           define lower bounds on the slack (control) variables
c
        do 120 i = 3,6
          y0(i) = zero
          y1(i) = zero
          do 110 j = -1,1
            ylb(j,i) = zero
 110      continue
 120    continue
c
C             Define objective function.
C
        iterm(1,1) = 0
        iterm(2,1) = 1
        iterm(3,1) = 0
        iterm(4,1) = -1
        coef(1) = 1.d0
c
        MAXMIN = -1
c
        init(1) = 1
c
      elseif(icl.eq.4) then
        init(1) = 2
      else
        print *,'invalid value for icl'
        stop
      endif
C
      IER = 0
C
      RETURN
      END
