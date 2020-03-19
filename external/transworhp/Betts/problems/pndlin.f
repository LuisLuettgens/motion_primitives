
      SUBROUTINE PNDLIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C             REF:  This is a mathematical pendulum problem posed by 
c                   Christof Buskens and Matthias Gerdts
C
C     ******************************************************************
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Arguments:
C
      INTEGER    IPHASE,NPHS,METHOD,NSTG(2),NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           CLBL(0:MXPCON)*80
      COMMON  /pndlCM/ ICLSS
c
      INCLUDE '../commons/odeprb.cmn'
      parameter (zero=0.d0)
C
C     ******************************************************************
C
C     ----PROBLEM TITLES AND LABELS (OPTIONAL)
C
      TITLE(1) = 'Pendulum Problem, (Buskens and Gerdts)'
c
      ngrid = itpgrd(iprob)
      method = itpmet(iprob)
      nstg(1) = itpstg(iprob)
      iclss = itpcls(iprob)
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),2;(la3),3;(la4),20')
cref      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
      call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      init = 1
c
      IF(ICLSS.EQ.1) THEN
c
        TITLE(2) = 'Index 3 DAE reduced to Index 1 DAE'
        STLBL(0) = 'TIME    Time'
        STLBL(1) = 'YSTATE1 State Variable 1'
        STLBL(2) = 'YSTATE2 State Variable 2'
        STLBL(3) = 'YSTATE3 State Variable 3'
        STLBL(4) = 'YSTATE4 State Variable 4'
        STLBL(5) = 'YSTATE5 Algebraic State Variable 5'
        STLBL(6) = 'UCNTRL  Control Variable'
C
C     ----NUMBER OF DIFFERENTIAL EQUATIONS (STATES)   
C
        NDE = 4
C
C     ----NUMBER OF ALGEBRAIC VARIABLES (CONTROLS)   
C
        NAV = 2
c
      elseIF(ICLSS.EQ.2) THEN
c
        TITLE(2) = 'Index 3 DAE reduced to an ODE'
        STLBL(0) = 'TIME    Time'
        STLBL(1) = 'YSTATE1 State Variable 1'
        STLBL(2) = 'YSTATE2 State Variable 2'
        STLBL(3) = 'YSTATE3 State Variable 3'
        STLBL(4) = 'YSTATE4 State Variable 4'
        STLBL(5) = 'YSTATE5 Algebraic State Variable 5'
        STLBL(6) = 'UCNTRL  Control Variable'
C
C     ----NUMBER OF DIFFERENTIAL EQUATIONS (STATES)   
C
        NDE = 5
C
C     ----NUMBER OF ALGEBRAIC VARIABLES (CONTROLS)   
C
        NAV = 1
c
      ENDIF
C
      NGRDGS = NGRID
      NCF(1) = NDE
c
c     ----default bounds
c
      do 120 j = -1,1
        do 110 i = 1,4
          ylb(j,i) = -5.d0
          yub(j,i) = 5.d0
 110    continue
        ylb(j,5) = -1.d0
        yub(j,5) = 15.d0
 120  continue
C
C     ----INITIAL TIME AND BOUNDARY CONDITION
C
      Y0(0) = 0.D0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
C
C     ----FINAL TIME AND BOUNDARY CONDITION
C
      TFINAL = 3.D0
      Y1(0) = TFINAL
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C     ----DEFINE INITIAL CONDITIONS FOR VARIABLES, Y.
C
      Y0(1) = 1.D0
      Y0(2) = 0.D0
      Y0(3) = 0.D0
      Y0(4) = 0.D0
      Y0(5) = 0.D0
      Y0(6) = 0.D0
C
C     ----DEFINE FINAL CONDITIONS FOR VARIABLES, Y.
C
      Y1(1) = 0.D0
      Y1(2) = 0.D0
      Y1(3) = 0.D0
      Y1(4) = 0.D0
      Y1(5) = 0.D0
      Y1(6) = 0.D0
C
C     ----FIX THE INITIAL STATES
C
      YLB(-1,1) = Y0(1)
      YUB(-1,1) = Y0(1)
      YLB(-1,2) = Y0(2)
      YUB(-1,2) = Y0(2)
      YLB(-1,3) = Y0(3)
      YUB(-1,3) = Y0(3)
      YLB(-1,4) = Y0(4)
      YUB(-1,4) = Y0(4)
      if(iclss.ne.1) then
        YLB(-1,5) = Y0(5)
        YUB(-1,5) = Y0(5)
      endif
C
C     ----FIX THE FINAL STATES
C
      YLB(1,1) = Y1(1)
      YUB(1,1) = Y1(1)
      YLB(1,3) = Y1(3)
      YUB(1,3) = Y1(3)
C
      NTERM =0
      NKON = 0
c
      if(iclss.eq.1) then
C
C     ----DEFINE path constraint
C
        ncf(2) = 1
        NTERM = NTERM + 1
        NKON = NKON + 1
        ITERM(1,NTERM) = NKON
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = -1
        COEF(NTERM) = 1.D0
        CLB(NKON) = zero
        CUB(NKON) = zero
        clBL(nkon) = 'PATHCN  Acceleration level Path constraint'
c
      endif
c
      IF(ICLSS.EQ.1) THEN
c
        ncf(3) = 1
C
C     ----DEFINE QUADRATURE AS THE OBJECTIVE TO BE MINIMIZED
C
        NTERM = NTERM + 1
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = -2
        COEF(NTERM) = 1.D0
C
        MAXMIN = -1
c
      elseIF(ICLSS.EQ.2) THEN
c
        ncf(3) = 1
C
C     ----DEFINE QUADRATURE AS THE OBJECTIVE TO BE MINIMIZED
C
        NTERM = NTERM + 1
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = -1
        COEF(NTERM) = 1.D0
C
        MAXMIN = -1
c
      ENDIF
C
      RETURN
      END
