
      SUBROUTINE QLININ(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                  INIT,MAXMIN,MXPARM,P0,PLB,PUB,plbl,
     &                  MXSTAT,Y0,Y1,YLB,YUB,STSKL,stlbl,MXPCON,CLB,CUB,
     &                  clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the quadratic-linear controller
C             problem.  
C
C             REF:  This is Chapter #5 in Bryson & Ho
C                   "APPLIED OPTIMAL CONTROL", P.148.
C
C             ICLASS = 1: Minimize final state.
C
C             ICLASS = 2: Minimize the control authority used to 
C                        drive a linear system of diff. eq. to zero.
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
      INCLUDE '../commons/odeprb.cmn'
C
      COMMON /QLCOM/ FMTRX(6,6),GMTRX(6,3),AMTRX(6,6),BMTRX(3,3),
     &               SFMTRX(6,6),ICL,NDEQ,NALG
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
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),2;(la3),3;(la4),20')
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
      call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      NPV = 0
      TITLE(1) = 'Quadratic-linear control problem'
C
      NAlg = 3
      NY  = 6
      IF (ICL.EQ.1) THEN
        NDEQ = NY
        TITLE(2) = 'Minimize quadrature of control'
      ELSEIF(ICL.EQ.2) THEN
        NDEQ = NY + 1
        TITLE(2) = 'Minimize integrated control'
      ELSEIF(ICL.EQ.3) THEN
        NDEQ = NY
        TITLE(2) = 'Minimize control quadrature with path constraint'
      ELSEIF(ICL.EQ.4) THEN
        NDEQ = 0
        ny = 0
        nalg = 1
        TITLE(2) = 'Minimize control quadrature--no dae'
      ENDIF
      NCF(1) = NDEQ
      NAV = nalg
c
      nterm =0
      nkon = 0
      npath = 0
C
C             Initial time and final time, fixed.
C
      Y0(0) = ZERO
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      if(icl.lt.4) then
        Y1(0) = ONEEP3
      else
        Y1(0) = ONE
      endif
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Load constant matrices describing diff. eq.
C
      DO J=1,NY
        DO I=1,NY
          AMTRX(I,J) = ZERO
          FMTRX(I,J) = ZERO
          SFMTRX(I,J) = ZERO
        enddo
      enddo
      FMTRX(1,4) = ONE
      FMTRX(2,5) = ONE
      FMTRX(3,6) = ONE
C
      DO J=1,NALG
        DO I=1,NY
          GMTRX(I,J) = ZERO
        enddo
      enddo
      GMTRX(4,1) = ONE
      GMTRX(5,2) = ONE
      GMTRX(6,3) = ONE
C
      DO J=1,NALG
        DO I=1,NALG
          BMTRX(I,J) = ZERO
        enddo
        BMTRX(J,J) = ONE
      enddo
C
C             Define initial conditions for state variables (X).
C
      ISGN = -1
      DO 75 I=1,NDEQ
        IF(I.LE.3) THEN
          Y0(I) = ONEEP3
        ELSEIF(I.LE.6) THEN
          Y0(I) = ONEEP1*ISGN
          ISGN = -ISGN
        ELSE
          Y0(I) = ZERO
        ENDIF
   75 CONTINUE
      DO 76 I=NDEQ+1,NDEQ+NALG
        Y0(I) = ZERO
   76 CONTINUE
C
C             Fix initial states.
C
      DO 90 I=1,NDEQ
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
   90 CONTINUE
C
C             NOTE: The control state and control bounds below
C                   are important for rapid convergence.
C
C             Bound integral of control state.
C
      IF (ICL.EQ.2) THEN
        YLB(0,7) = ZERO
        YUB(0,7) = 100.D0
      ENDIF
C
C             Bound accelerations (controls).
C
      IF(ICL.LT.4) THEN
        CTLBND = ONE
      ELSE
        CTLBND = 2.D0
      ENDIF
      DO 93 I=NDEQ+1,NDEQ+NALG
        YLB(0,I) = -CTLBND
        YUB(0,I) = CTLBND
   93 CONTINUE
C
C             Define terminal conditions for state variables (X).
C
      DO 95 I=1,NDEQ+NALG
        Y1(I) = ZERO
   95 CONTINUE
C
C             Fix final states.
C
      DO 110 I=1,NY
        YLB(1,I) = Y1(I)
        YUB(1,I) = Y1(I)
  110 CONTINUE
c
      if(icl.eq.3) then
        ncf(2) = ncf(2) + 1
        npath = npath - 1
        nterm = nterm + 1
        nkon = nkon + 1
        clb(nkon) = -1.d4
        cub(nkon) = 1.d4
        iterm(1,nterm) = nkon
        iterm(2,nterm) = 1
        iterm(3,nterm) = 0
        iterm(4,nterm) = npath
        coef(nterm) = one
      endif
c
C             Define objective function.
C
      nterm = nterm + 1
      iterm(1,nterm) = 0
      iterm(2,nterm) = 1
      if(icl.eq.2) then
        ITERM(3,1) = 1
        ITERM(4,1) = NDEQ
      else
        ncf(3) = ncf(3) + 1
        npath = npath - 1
        iterm(3,nterm) = 0
        iterm(4,nterm) = npath
      endif
      coef(nterm) = one
C
      MAXMIN = -1
C
      INIT = 1
C
C             Load state scale weights for output.
C
      DO 150 I=0,NDEQ+NALG
        STSKL(I,1) = ONE
  150 CONTINUE
C
      IER = 0
C
      RETURN
      END
