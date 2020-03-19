      SUBROUTINE ARAOIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C             REF:  This is a problem suggested by Anil V. Rao
C
C     ******************************************************************
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Arguments:
C
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT(2),MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           CLBL(0:MXPCON)*80
      COMMON  /ARAOCM/ TFINAL,ICLSS,NGRDGS
c
      INCLUDE '../commons/odeprb.cmn'
C
C     ******************************************************************
C
C     ----PROBLEM TITLES AND LABELS (OPTIONAL)
C
      TITLE(1) = 'A. V. Rao Test Problem'
C
      ICLSS = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
      IF(TPSEED(IPROB).GE.0.D0) THEN
        INIT(1) = 1
      ELSE
        INIT(1) = 2
      ENDIF
c
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),20')
      call STMSEQ('(trp),2;(la3),20')
C
C     ----NUMBER OF ALGEBRAIC VARIABLES (CONTROLS)   
C
      NAV = 1
C
      IF(ICLSS.EQ.1) THEN
c
        TITLE(2) = 'Lagrange Formulation'
        STLBL(0) = 'TIME    Time'
        STLBL(1) = 'XSTATE  State Variable'
        STLBL(2) = 'UCNTRL  Control Variable'
C
C     ----NUMBER OF DIFFERENTIAL EQUATIONS (STATES)   
C
        NDE = 1
C
      ELSEIF(ICLSS.EQ.2) THEN
c
        TITLE(2) = 'Mayer Formulation'
        STLBL(0) = 'TIME    Time'
        STLBL(1) = 'XSTAT1  State Variable 1'
        STLBL(2) = 'XSTAT2  State Variable 2'
        STLBL(3) = 'UCNTRL  Control Variable'
C
C     ----NUMBER OF DIFFERENTIAL EQUATIONS (STATES)   
C
        NDE = 2
C
      ELSEIF(ICLSS.EQ.3) THEN
c
        TITLE(2) = 'Indirect Formulation'
        STLBL(0) = 'TIME    Time'
        STLBL(1) = 'XSTAT1  State Variable 1'
        STLBL(2) = 'ADJ1    Adjoint Variable 1'
C
C     ----NUMBER OF DIFFERENTIAL EQUATIONS 
C
        NDE = 2
        NAV = 0
C
      ENDIF
C
      NGRDGS = NGRID
      NCF(1) = NDE
C
C     ----INITIAL TIME AND BOUNDARY CONDITION
C
      Y0(0) = 0.D0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
C
C     ----FINAL TIME AND BOUNDARY CONDITION
C
      TFINAL = 10000.D0
      Y1(0) = TFINAL
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C     ----DEFINE INITIAL CONDITIONS FOR STATE VARIABLES, Y.
C
      Y0(1) = 1.D0
C
C     ----DEFINE FINAL CONDITIONS FOR STATE VARIABLES, Y.
C
      Y1(1) = 1.5D0
C
C     ----FIX THE INITIAL AND FINAL STATES
C
      YLB(-1,1) = Y0(1)
      YUB(-1,1) = Y0(1)
      YLB(1,1) = Y1(1)
      YUB(1,1) = Y1(1)
C
C     ----ADDITIONAL STATE FOR CLASS = 2
C
      IF(ICLSS.EQ.2) THEN
        Y0(2) = 0.D0
        Y1(2) = 1.0D0
        YLB(-1,2) = Y0(2)
        YUB(-1,2) = Y0(2)
c       state scaling dramatically improves constraint solve (11/6/04)
        stskl(2,2) = 1.d-2
      ENDIF
C
      IF(ICLSS.LE.2) THEN
C
C     ----DEFINE INITIAL AND FINAL CONTROL ANGLES.
C
        Y0(NDE+1) = 0.D0
        Y1(NDE+1) = 0.D0
C
      ELSE
C
C     ----DEFINE INITIAL AND FINAL ADJOINT
C
        Y0(2) = 8.2946869E-01
        Y1(2) = -1.4028584E+01
C
      ENDIF
C
      NTERM =0
      IF(ICLSS.EQ.1) THEN
C
C     ----DEFINE QUADRATURE AS THE OBJECTIVE TO BE MINIMIZED
C
        NCF(3) = 1
        NTERM = NTERM + 1
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = -1
        COEF(NTERM) = 1.D0
C
        MAXMIN = -1
      ELSEIF(ICLSS.EQ.2) THEN
C
C     ----DEFINE SECOND STATE AS THE OBJECTIVE TO BE MINIMIZED
C
        NTERM = NTERM + 1
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 1
        ITERM(4,NTERM) = NDE
        COEF(NTERM) = 1.D0
C
        MAXMIN = -1
      ELSEIF(ICLSS.EQ.3) THEN
        MAXMIN = 0
      ENDIF
C
      RETURN
      END
