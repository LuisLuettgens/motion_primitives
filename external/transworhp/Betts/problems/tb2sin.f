
      SUBROUTINE TB2SIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C             "Optimal Control of Treatments in a Two-Strain Tuberculosis Model"
C             E. Jung, S. Lenhart, Z. Feng,  Discrete and Continuous Dynamical
c             Systems--Series B.,  Volume 2, Number 4, November 2002, pp 473-482
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
c
      parameter  (beta1 = 13.d0, beta2 = 13.d0)
      parameter  (xmu = .0143d0, done  = 0.d0 , dtwo   = 0.d0  )
      parameter  (xk1 = .5d0,    xk2   = 1.d0 , rone   = 2.d0  )
      parameter  (rtwo = 1.d0,   smlp  = .4d0 , smlq   = .1d0  )
      parameter  (capN = 30.d3,  bstar  = .029d0 )
      parameter  (Bone  = 50.d0, Btwo   = 500.d0)
      parameter  (capLam = xmu*capN)
      parameter  (Szero  = 76.d0*capN/120.d0)
      parameter  (capL10 = 36.d0*capN/120.d0)
      parameter  (capI10 =  4.d0*capN/120.d0)
      parameter  (capL20 =  2.d0*capN/120.d0)
      parameter  (capI20 =  1.d0*capN/120.d0)
      parameter  (capT0  =  1.d0*capN/120.d0)
      parameter  (zero = 0.d0, one = 1.d0)
c
      INCLUDE '../commons/odeprb.cmn'
c
C     ==================================================================
c
      TITLE(1) = 'Optimal Treatment in a Two-Strain Tuberculosis Model'
      TITLE(2) = 'Jung, Lenhart, and Feng'
c
      ICLS = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
      seedtp = tpseed(iprob)
c
cold      call STMSEQ('(la2),2;(la3),20')
cref      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
      call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
c
      init(1) = 1
c
      nphs = 1
      nterm = 0
      nkon = 0
      bigbnd = one/hdmcon(5)
      nde = 0
      nav = 0
c
C     ==================================================================
C     ========Differential Variables====================================
C     ==================================================================
c
C             Set dynamic variable labels.
c
      STLBL(0) = 'TIME    Time (years)'
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'capS    Susceptible'
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'capL1   Latent, Infected, Typ. TB, not infectious'
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'capI1   Infectious, Typical TB'
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'capL2   Latent, Infected, Res.  TB, not infectious'
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'capI2   Infectious, Resistant TB'
c
c     ------------------------------
c
      nde = nde + 1
      STLBL(nde) = 'capT    Treated (effectively)'
c  
      NCF(1) = NDE
c
c     ==================================================================
C     ========Algebraic Variables=======================================
C     ==================================================================
c
      bndlwr = +.05d0
      bndupr = +.95d0
c
      nav = nav + 1
      STLBL(nde+nav) = 'UONE    Case Finding Control'
c
      call dfill(3,bndlwr,ylb(-1,nde+nav),1)
      call dfill(3,bndupr,yub(-1,nde+nav),1)
c
c     ------------------------------
c
      nav = nav + 1
      STLBL(nde+nav) = 'UTWO    Case Holding Control'
c
      call dfill(3,bndlwr,ylb(-1,nde+nav),1)
      call dfill(3,bndupr,yub(-1,nde+nav),1)
c
c     ==================================================================
C     ========Boundary Conditions=======================================
C     ==================================================================
c
c         load initial guess 
c
      y0(0) =  zero
      y0(1) =  Szero 
      y0(2) =  capL10
      y0(3) =  capI10
      y0(4) =  capL20
      y0(5) =  capI20
      y0(6) =  capT0 
      y0(7) =  zero
      y0(8) =  zero
c
      y1(0) =  5.d0
      call dcopy(nde+nav,y0(1),1,y1(1),1)
c
      do kk = 0,nde
        ylb(-1,kk) = y0(kk)
        yub(-1,kk) = y0(kk)
      enddo
c
      ylb(+1,0) = y1(0)
      yub(+1,0) = y1(0)
c
c     ==================================================================
C     =========Objective Function=======================================
C     ==================================================================
C
c
      MAXMIN = -1
c
      CLBL(0) = 'TBGOAL  TB Treatment Goals'
C
C           OBJECTIVE FUNCTION DEFINITION.
C
      ncf(3) = ncf(3) + 1
      NTERM = NTERM + 1
      ITERM(1,NTERM) = 0
      ITERM(2,NTERM) = IPHASE
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = -ncf(2) -ncf(3)
      COEF(NTERM) = one
C
      RETURN
      END
