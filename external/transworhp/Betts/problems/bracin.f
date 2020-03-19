

      SUBROUTINE BRACIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             INITIALIZE BRACHISTOCHRONE CONTROL EXAMPLE.
C
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,ONEEP1=1.D+01,ONEEP3=1.D+03)
C
      include '../commons/odeprb.cmn'
C
C  ARGUMENTS:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      character  title(3)*60,plbl(mxparm)*80,stlbl(0:mxstat)*80,
     &           clbl(0:mxpcon)*80
      COMMON /BRACCM/ HBRAC,ICLBRC,NGRDBR
C
C  LOCAL:
C
C
C     ******************************************************************
C
C         MAXIMUM NUMBER OF PHASES
C
      NPHS = 1
C
C         PROBLEM CLASS
C
      ICLBRC = ITPCLS(IPROB)
C
C         NUMBER OF GRID POINTS
C
C
      NGRID = ITPGRD(IPROB)
      NGRDBR = NGRID
C
C         DISCRETIZATION METHOD (TRAPEZOIDAL=2)
C
      METHOD = ITPMET(IPROB)
C
C         NUMBER OF STAGES---NOT USED FOR TRAPEZOIDAL METHOD
C
      NSTG = ITPSTG(IPROB)
C
C         PROBLEM SEED VALUE--SELECT INITIALIZATION AND PROBLEM TYPE
C
      SEEDBR = TPSEED(IPROB)
      HBRAC = ABS(SEEDBR)
c
      select case (iclbrc)
      case(1)
cold        call STMSEQ('(la2),2;(la3),3;(la4),20')
cref        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
        call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      case(2)
cref        call STMSEQ('(la2),2;(la3),20')
        call STMSEQ('(trp),2;(la3),20')
      end select
C
      IF(SEEDBR.LT.ZERO) THEN
C
C         (INIT = 2 --- USE SUBROUTINE BRACIG)
C
        INIT = 2
        T1 = SQRT(HDMCON(12)/GRAV)
C
      ELSE
C
C         INITIAL GUESS (INIT=1 --- LINEAR) 
C
        INIT = 1
        T1 = 2.D0
C
      ENDIF
C
C         PROBLEM TITLES
C
      TITLE(1) = 'BRACHISTOCHRONE CONTROL EXAMPLE'
      IF(ICLBRC.EQ.1) THEN
        TITLE(2) = 'UNCONSTRAINED ANALYTIC SOLUTION'
      ELSEIF(ICLBRC.EQ.2) THEN
        TITLE(2) = 'STATE VARIABLE INEQUALITY CONSTRAINT'
      ELSEIF(ICLBRC.EQ.3) THEN
        TITLE(2) = 'STATE VARIABLE INEQUALITY--REDUCED INDEX'
      ELSE
        PRINT *,'INCORRECT VALUE FOR ICLBRC =',ICLBRC
        STOP
      ENDIF
C
C         NUMBER OF ALGEBRAIC VARIABLES (# CONTROLS)
C
      NAV = 1
C
C         NUMBER OF PARAMETRIC VARIABLES
C
      NPV = 0
C
C         NUMBER OF DIFFERENTIAL EQUATIONS
C
      NDE = 3
      NCF(1) = NDE
C
C         NUMBER OF ALGEBRAIC EQUATIONS
C
      IF(ICLBRC.GE.2) THEN
        NAE = 1
        NCF(2) = NAE
      ENDIF
C
C         FIXED INITIAL TIME AND VALUE
C
      T0 = ZERO
      Y0(0) = T0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
C
C         FREE FINAL TIME AND GUESS
C
      Y1(0) = T1
C
C         LOWER BOUND FOR FINAL TIME
C
      YLB(1,0) = ZERO
C
C         DEFINE INITIAL CONDITIONS FOR DYNAMIC (STATE+CONTROL) VARIABLES 
C
      GRAV  = 32.174D0
      CR2D  = HDMCON(16)
      Y0(1) = ZERO
      Y0(2) = ZERO
      Y0(3) = ZERO
      Y0(4) = ZERO
C
C         FIX INITIAL STATES.
C
      YLB(-1,1) = Y0(1)
      YUB(-1,1) = Y0(1)
      YLB(-1,2) = Y0(2)
      YUB(-1,2) = Y0(2)
      YLB(-1,3) = Y0(3)
      YUB(-1,3) = Y0(3)
C
C         BOUND STATES 
C
      YLB(0,1) = ZERO
      YUB(0,1) = ONEEP1
      YLB(0,2) = ZERO
      YUB(0,2) = ONEEP1
      YLB(0,3) = ZERO
      YUB(0,3) = ONEEP1
c
      BIGBND = ONE/HDMCON(5)
C
C         BOUND CONTROL (VAR #4) AT BEGINNING (=-1), DURING (=0), AND
C         END OF PHASE (=1)
C
      YLB(-1,4) = ZERO
      YUB(-1,4) = (90.D0/CR2D)
      YLB(0,4) = ZERO
      YUB(0,4) = (90.D0/CR2D)
      YLB(1,4) = ZERO
      YUB(1,4) = (90.D0/CR2D)
C
C         DEFINE TERMINAL CONDITIONS FOR STATE VARIABLES (X).
C
      Y1(1) = ONE
      Y1(2) = ONE
      Y1(3) = ONE
      Y1(4) = ZERO
C
C         FIX FINAL STATE (HORIZONTAL DISTANCE)
C
      YLB(1,1) = Y1(1)
      YUB(1,1) = Y1(1)
C
      NTERM = 0
      NKON = 0
C
C         STATE VARIABLE PATH CONSTRAINT
C
      IF(ICLBRC.GE.2) THEN
        NTERM = NTERM + 1
        NKON = NKON + 1
        ITERM(1,NTERM) = NKON
        ITERM(2,NTERM) = IPHASE
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = -1
        COEF(NTERM) = ONE
        CLB(NKON) = -BIGBND
        CUB(NKON) = ZERO
      ENDIF
C
C         MINIMIZE THE FINAL TIME
C
      MAXMIN = -1
C
C         DEFINE OBJECTIVE.
C
      NTERM = NTERM + 1
      ITERM(1,NTERM) = 0
      ITERM(2,NTERM) = 1
      ITERM(3,NTERM) = 1
      ITERM(4,NTERM) = 0
      COEF(NTERM) = ONE
C
C             LOAD STATE SCALE WEIGHTS FOR OUTPUT.
C
      DO 150 I=1,1+NDE+NAV
        STSKL(I,1) = ONE
  150 CONTINUE
C
      IER = 0
C
      RETURN
      END
