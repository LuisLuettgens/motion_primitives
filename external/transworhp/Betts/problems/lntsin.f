



      SUBROUTINE LNTSIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the linear tangent steering collocation
C             problem for the state/adjoint or the state/control
C             formulation.
C
C             REF:  This is Problem #9 in Bryson & Ho
C                   "APPLIED OPTIMAL CONTROL", P.82-83.
C                   The analytic answer is listed in the comments.
C
C             ICLASS = 1: Solve optimal control problem by solving
C                        nonlinear system of equations using adjoints.
C
C             ICLASS = 2: Solve optimal control problem by minimization
C                        using controls directly as variables.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FIVE=5.D0)
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
C
C  Local:
      COMMON /LNTSCM/ ATHRUS,CPI2,ICL
C
      INCLUDE '../commons/odeprb.cmn'
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
      TITLE(1) = 'Linear Tangent Steering ascent'
C
      NPV = 0
      INIT(1) = 1
      STLBL(0) = 'TIME    Time'
      STLBL(1) = 'RANGE   Range'
      STLBL(2) = 'ALTITUDE Altitude'
      STLBL(3) = 'RNG-RATE Range Rate'
      STLBL(4) = 'ALT-RATE Altitude Rate'
      PLBL(1) = 'TFINAL   Final time'
      IF (ICL.EQ.1) THEN
          TITLE(2) = 'Solve by adjoints'
          NDE = 8
          NAV = 0
          STLBL(1) = 'RANGE   Range'
          STLBL(2) = 'ALTITUDE Altitude'
          STLBL(3) = 'RNG-RATE Range Rate'
          STLBL(4) = 'ALT-RATE Altitude Rate'
          STLBL(5) = 'ADJ-RNG Range Adjoint Variable'
          STLBL(6) = 'ADJ-ALT Altitude Adjoint Variable'
          STLBL(7) = 'ADJ-RRAT Range Rate Adjoint Variable'
          STLBL(8) = 'ADJ-ARAT Altitude Rate Adjoint Variable'
c
cold          call STMSEQ('(la2),2;(la3),20')
cold          call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref          call STMSEQ('(la2),2;(la3),3;(la4),20')
          call STMSEQ('(trp),2;(la3),3;(la4),20')
      ELSEIF (ICL.EQ.2) THEN
          TITLE(2) = 'Solve by controls'
          NDE = 4
          NAV = 1
          STLBL(5) = 'PITCH   Pitch angle (deg)'
c
cold          call STMSEQ('(la2),2;(la3),20')
cold          call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref          call STMSEQ('(la2),2;(la3),3;(la4),20')
          call STMSEQ('(trp),2;(la3),3;(la4),20')
      elseif (icl.eq.3) then
          TITLE(2) = 'Right Hand Side Discontinuity'
          NDE = 4
          NAV = 1
          STLBL(5) = 'PITCH   Pitch angle (deg)'
c
      ELSEIF (ICL.EQ.4) THEN
          INIT(1) = 6
          INIT(2) = 2
          TITLE(2) = 'Explicit Linear Tangent Parameterization'
          NDE = 4
          NAV = 0
          STLBL(5) = 'PITCH   Pitch angle (deg)'
      ELSE
          IER = -11
          RETURN
      ENDIF
      NCF(1) = NDE
C
C             Initial time (fixed) and guess for final time.
C
      Y0(0) = ZERO
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      T1 = .554570878337D0
      Y1(0) = T1
      YLB(1,0) = ZERO
C
C             Define initial conditions for state variables, Y.
C                 Y(1) = Range.
C                 Y(2) = Altitude.
C                 Y(3) = Horizontal velocity.
C                 Y(4) = Vertical velocity.
C
C             Declare initial values as fixed.
C
      CR2D = HDMCON(16)
      CPI2 = HDMCON(12)/TWO
      ATHRUS = 100.D0
      Y0(1) = ZERO
      Y0(2) = ZERO
      Y0(3) = ZERO
      Y0(4) = ZERO
      DO 20 I=1,4
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
   20 CONTINUE
C
C             Define constants for closed form solution.
C
      B0 = 54.6263551908D0/CR2D
      SINB0 = SIN(B0)
      COSB0 = COS(B0)
      TANB0 = SINB0/COSB0
      SECB0 = ONE/COSB0
      C = TWO*TANB0/T1
C
C          T = T1
C          SECB0 = 1./COS(B0)
C
C          TANB = TANB0 - C*T
C          BETA = ATAN(TANB)
C          SECB = 1./COS(BETA)
C          ARG1 = ALOG((TANB0+SECB0)/(TANB+SECB))
C          ARG2 = SECB0 - SECB
C          ARG3 = TANB0 - TANB
C          Y(1) = ATHRUS*(ARG2 -TANB*ARG1)/C**2
C          Y(2) = ATHRUS*(ARG3*SECB0-ARG2*TANB-ARG1)/(2.*C**2)
C          Y(3) = ATHRUS*ARG1/C
C          Y(4) = ATHRUS*ARG2/C
C          Y(5) = 0.
C          Y(6) = -2.*SIN(B0)/(ATHRUS*T1)
C          Y(7) = -COS(B0)/ATHRUS
C          Y(8) = -SIN(B0)*(1.-2.*T/T1)/ATHRUS
C
C
C             Define terminal conditions for state variables.
C
      TANB = TANB0 - C*T1
      BETA = ATAN(TANB)
      SECB = ONE/COS(BETA)
      ARG1 = LOG( (TANB0+SECB0)/(TANB+SECB) )
      Y1(1) = 12.4778447625D0
      Y1(2) = FIVE
      Y1(3) = ATHRUS*ARG1/C
      Y1(4) = ZERO
C
      IF (ICL.EQ.1) THEN
C
C             Define initial adjoints.
C
          Y0(5) = ZERO
          Y0(6) = -TWO*SINB0/(ATHRUS*T1)
          Y0(7) = -COSB0/ATHRUS
          Y0(8) = -SINB0/ATHRUS
C
C             Define final adjoints as fixed.
C
          Y1(5) = ZERO
          Y1(6) = Y0(6)
          Y1(7) = Y0(7)
          Y1(8) = -Y0(8)
c
c           fix final states and first adjoint
c
          ylb(+1,0) = .01d0
c
          ylb(1,2:5) = y1(2:5)
          yub(1,2:5) = y1(2:5)
C
C             Define optimality constraint.
C
          ITERM(1,1) = 1
          ITERM(2,1) = 1
          ITERM(3,1) = 1
          ITERM(4,1) = -1
          COEF(1) = ONE
          CLB(1) = ZERO
          CUB(1) = ZERO
C
          NPF(2) = 1
          MAXMIN = 0
      elseif(icl.eq.4) then
c
c     ==================================================================
C     ========Parametric Variables======================================
C     ==================================================================
c
        npv = npv + 1
        plbl(npv) = 'PONE    Linear Tangent Parameter 1'
        plb(npv) = 0.d0
        pub(npv) = 10.d0
        p0(npv) = 1.d0
c
        npv = npv + 1
        plbl(npv) = 'PTWO    Linear Tangent Parameter 2'
        plb(npv) = 0.d0
        pub(npv) = 10.d0
        p0(npv) = 1.d0
c
c           fix final states and first adjoint
c
        ylb(+1,0) = .01d0
c
        ylb(1,2:4) = y1(2:4)
        yub(1,2:4) = y1(2:4)
c
        MAXMIN = 0
c
      ELSE
C
C             Define initial and final control angles.
C
          Y0(5) = B0
          BF = ATAN(TANB0-C*T1)
          Y1(5) = BF
C
C             Fix final states 2-4.
C
          DO 70 I=2,4
            YLB(1,I) = Y1(I)
            YUB(1,I) = Y1(I)
   70     CONTINUE
C
C             Bound control magnitude on phase by 90.
C
          DO 75 I=-1,1
            YLB(I,5) = -CPI2
            YUB(I,5) =  CPI2
   75     CONTINUE
C
C             Define final time as objective.
C
          ITERM(1,1) = 0
          ITERM(2,1) = 1
          ITERM(3,1) = 1
          ITERM(4,1) = 0
          COEF(1) = ONE
C
          MAXMIN = -1
      ENDIF
C
C             Load state scale weights for output.
C
      DO 150 I=0,NDE
        STSKL(I,1) = ONE
  150 CONTINUE
      DO 152 I=NDE+1,NDE+NAV
        STSKL(I,1) = CR2D
  152 CONTINUE
C
      RETURN
      END
