
      SUBROUTINE HEATIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the distributed parameter optimal control
C             of a heat diffusion process. References:
c
c             icl = 1  Betts and Citron, "Approximate Optimal Control
c                      of Distributed Parameter Systems," AIAA Journal, Vol 10,
c                      No. 1, jan 1972, pp 19-23
c
c             icl = 2  Matthias Heinkenschloss, "Projected Sequential
c                      Quadratic Programming Methods," SIAM J. Optimization, 
c                      Vol 6, No. 2, pp 373-417, May, 1996.
c
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,C5P0=5.D0)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT),CLB(MXPCON),CUB(MXPCON),COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           CLBL(0:MXPCON)*80
C
      COMMON /HEATC1/ TFINAL,VUPR,CAPH,VLWR,QSUBA,GAMMA,HSTEP,
     &             PI,ACOEF1,ACOEF2,ACOEF3,ACOEF4,GCOEF,RHO,GAMRLX
      COMMON /HEATC2/ ICL,NXGRID
C
      INCLUDE '../commons/odeprb.cmn'
C
      CHARACTER*60 LABEL
      DATA LABEL(1:60) / ' '/
C
C     ******************************************************************
C
      ICL = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
c
cok      call STMSEQ('(la2),2;(la3),20')
cref      call STMSEQ('(la2),-2;(la3),-3;(la4),-4;(la5),-20')
      call STMSEQ('(trp),-2;(la3),-3;(la4),-4;(la5),-20')
C
C         DEFINE THE SPATIAL DISCRETIZATION AND VARIOUS PROBLEM 
C         PARAMETERS
C
      NXGRID = TPSEED(IPROB)
      bigbnd = one/hdmcon(5)
c
      if(icl.eq.1) then
        TFINAL = .2D0
        VUPR = ONE
        CAPH = 10.D0
        VLWR = ZERO
        QSUBA = .2D0
        GAMMA = .04D0
        lcntrl = nxgrid + 1
      elseif(icl.eq.2) then
        acoef1 = 4.0d0
        acoef2 = one
        acoef3 = 4.0d0
        acoef4 = -one
        gcoef  = one
        rho    = -one
        tfinal = .5d0
        gamrlx = 1.0d-3
        pi     = hdmcon(12)
        vlwr   = -bigbnd
        vupr   = .1d0
        lcntrl = nxgrid
      else
        print *,'icl is wrong'
        stop
      endif
C         ----SPATIAL DISCRETIZATION STEPSIZE
      HSTEP = ONE/REAL(NXGRID-1)
C
C             OBTAIN NUMBER OF GRID POINTS FOR PROBLEM.
C
      NPHS = 1
      if(icl.eq.1) then
        TITLE(1) = 'Heat Diffusion Process (Betts)'
      else
        TITLE(1) = 'Heat Diffusion Process (Heinkenschloss)'
      endif
      TITLE(2) = 'PDE Solved Using Method of Lines'
      STLBL(0) = 'TIME    Time'
C
      DO 108 I = 1,NXGRID
        LABEL(1:43) = 'Q       Temperature at x =                 '
        WRITE(LABEL(2:3),'(I2)') I
        label(1:3) = adjustl(label(1:3))
        XLOCI = REAL(I-1)/REAL(NXGRID-1)
        WRITE(LABEL(27:43),'(G16.6)') XLOCI
        STLBL(I) = LABEL
 108  CONTINUE
c
      if(icl.eq.1) then
        STLBL(NXGRID+1) = 'W(T)    Input temperature w(t)'
        STLBL(NXGRID+2) = 'V(T)    Fuel flow v(t)'
        LABEL(1:43) = 'Q       Temperature at x =                 '
        WRITE(LABEL(2:3),'(I2)') 00
        label(1:3) = adjustl(label(1:3))
        XLOCI = REAL(-1)/REAL(NXGRID-1)
        WRITE(LABEL(27:43),'(G16.6)') XLOCI
        STLBL(NXGRID+3) = LABEL
        LABEL(1:43) = 'Q       Temperature at x =                 '
        WRITE(LABEL(2:3),'(I2)') NXGRID+1
        label(1:3) = adjustl(label(1:3))
        XLOCI = REAL(NXGRID+1-1)/REAL(NXGRID-1)
        WRITE(LABEL(27:43),'(G16.6)') XLOCI
        STLBL(NXGRID+4) = LABEL
C
        NPV = 0
        INIT = 1
        NDE = NXGRID + 1
        NAV = 3
        NCF(1) = NDE
      else
        STLBL(NXGRID+1) = 'U(T)    Input temperature u(t)'
        LABEL(1:43) = 'Q       Temperature at x =                 '
        WRITE(LABEL(2:3),'(I2)') 00
        label(1:3) = adjustl(label(1:3))
        XLOCI = REAL(-1)/REAL(NXGRID-1)
        WRITE(LABEL(27:43),'(G16.6)') XLOCI
        STLBL(NXGRID+2) = LABEL
        LABEL(1:43) = 'Q       Temperature at x =                 '
        WRITE(LABEL(2:3),'(I2)') NXGRID+1
        label(1:3) = adjustl(label(1:3))
        XLOCI = REAL(NXGRID+1-1)/REAL(NXGRID-1)
        WRITE(LABEL(27:43),'(G16.6)') XLOCI
        STLBL(NXGRID+3) = LABEL
C
        NPV = 0
        INIT = 1
        NDE = NXGRID 
        NAV = 3
        NCF(1) = NDE
      endif
C
C         INITIALIZE THE CONSTRAINT DEFINITION PARAMETERS
C
      NCF(2) = 0
      NTERM =0
      NKON = 0
      NPATH = 0
C
C
C         INITIAL AND FINAL TIME FOR THE PROCESS
C
      Y0(0) = ZERO
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      Y1(0) = TFINAL
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C         DEFINE INITIAL CONDITIONS FOR STATE VARIABLES, Y.
C
      if(icl.eq.1) then
        CALL DFILL(NDE,ZERO,Y0(1),1)
      else
        do 109 i = 1,nxgrid
          xsubi = real(i-1)/real(nxgrid-1)
          y0(i) = two + cos(pi*xsubi)
 109    continue
      endif
C
C         DEFINE INITIAL CONDITIONS FOR CONTROL VARIABLES.
C
      Y0(NDE+1) = zero
      Y0(NDE+2) = Y0(1)
      Y0(NDE+3) = Y0(NXGRID)
C
C         FIX THE INITIAL VALUES FOR THE TEMPERATURE STATES
C
      DO 110 I=1,NDE
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
  110 CONTINUE
C
C         DEFINE TERMINAL CONDITIONS FOR STATE VARIABLES.
C
      if(icl.eq.1) then
        CALL DFILL(NDE,QSUBA,Y1(1),1)
      else
        CALL DFILL(NDE,zero,Y1(1),1)
      endif
C
C         DEFINE TERMINAL CONDITIONS FOR CONTROL VARIABLES.
C
      if(icl.eq.1) then
        Y1(NDE+1) = Y1(NDE)
        Y1(NDE+2) = Y1(1)
        Y1(NDE+3) = Y1(NXGRID)
      else
        Y1(NDE+1) = zero
        Y1(NDE+2) = Y1(1)
        Y1(NDE+3) = Y1(NXGRID)
      endif
C
C         BOUNDS ON THE CONTROL VARIABLE
C
      DO 120 I = -1,1
        YLB(I,lcntrl+1) = VLWR
        YUB(I,lcntrl+1) = VUPR
 120  CONTINUE
C
C         PATH CONSTRAINTS 
C
C         ------
C
C         BOUNDARY CONDITION AT X = 0
C
      NCF(2) = NCF(2) + 1
      NPATH = NPATH - 1
      NTERM = NTERM + 1
      NKON = NKON + 1
      CLBL(NKON) = 'DQDX0   Boundary condition on dq/dx at x = 0'
      CLB(NKON) = ZERO
      CUB(NKON) = ZERO
      ITERM(1,NTERM) = NKON
      ITERM(2,NTERM) = 1
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = NPATH
      COEF(NTERM) = ONE
C
C         BOUNDARY CONDITION AT X = 1
C
      NCF(2) = NCF(2) + 1
      NPATH = NPATH - 1
      NTERM = NTERM + 1
      NKON = NKON + 1
      CLBL(NKON) = 'DQDX1   Boundary condition on dq/dx at x = 1'
      CLB(NKON) = ZERO
      CUB(NKON) = ZERO
      ITERM(1,NTERM) = NKON
      ITERM(2,NTERM) = 1
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = NPATH
      COEF(NTERM) = ONE
C
      if(icl.eq.1) then
C
CONEPNT --->     FORMULATION USING ONE POINT FUNCTION 
CONEPNTC         DEFINE OBJECTIVE COMPUTED BY POINT FUNCTION.
CONEPNTC
CONEPNT      NPF(2) = NPF(2) + 1
CONEPNT      NTERM = NTERM + 1
CONEPNT      CLBL(0) = 'TMPDEV  Temperature deviation integral'
CONEPNT      ITERM(1,NTERM) = 0
CONEPNT      ITERM(2,NTERM) = 1
CONEPNT      ITERM(3,NTERM) = +1
CONEPNT      ITERM(4,NTERM) = -1
CONEPNT      COEF(1) = ONE
C         DEFINE OBJECTIVE COMPUTED BY NXGRID POINT FUNCTIONS.
C
        DO 130 I = 1,NXGRID
          NPF(2) = NPF(2) + 1
          NTERM = NTERM + 1
          ITERM(1,NTERM) = 0
          ITERM(2,NTERM) = 1
          ITERM(3,NTERM) = +1
          ITERM(4,NTERM) = -I
          COEF(1) = ONE
 130    CONTINUE
        CLBL(0) = 'TMPDEV  Temperature deviation integral'
c
      else
C
C         Temperature deviation integral
C
        NCF(3) = NCF(3) + 1
        NPATH = NPATH - 1
        NTERM = NTERM + 1
        CLBL(0) = 'TMPDEV  Temperature deviation integral'
        ITERM(1,NTERM) = 0
        ITERM(2,NTERM) = 1
        ITERM(3,NTERM) = 0
        ITERM(4,NTERM) = NPATH
c
      endif
C
      MAXMIN = -1
C
      RETURN
      END
