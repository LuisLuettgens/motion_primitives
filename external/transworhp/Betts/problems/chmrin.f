


      SUBROUTINE CHMRIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 clbl,MXTERM,COEF,ITERM,TITLE,IER)
C
C             Initialize the chemical reactor problem
C
C             REF:  pp. 123-131, Elements of Optimal Control, Holt,
c                   Rinehart and Wintson, Inc., 1969, Stephen J. Citron
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,FIVE=5.D0)
C
C  Arguments:
      INTEGER    IPHASE,NPHS,METHOD,NSTG(2),NCF(3),NPF(2),NPV,NAV,NGRID,
     &           INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,
     &           ITERM(4,MXTERM),IER
      DIMENSION  P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     &           Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     &           STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     &           COEF(MXTERM)
      CHARACTER  TITLE(3)*60,PLBL(MXPARM)*80,STLBL(0:MXSTAT)*80,
     &           clbl(0:mxpcon)*80
C
C  Local:
      COMMON /CHMRCM/ alwr,aupr,tfinal,xzero,yzero,rhocof,expk
C
      INCLUDE '../commons/odeprb.cmn'
C
C     ******************************************************************
C
C             Obtain number of grid points for problem.
C
      NPHS = 1
      ICASE = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG(1) = ITPSTG(IPROB)
c
      select case(icase)
      case(1)
cold        call STMSEQ('(la2),2;(la3),3;(la4),20')
cref        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
        call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      case(2)
cold        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref        call STMSEQ('(la2),2;(la3),3;(la4),20')
        call STMSEQ('(trp),2;(la3),3;(la4),20')
      case(3)
cold        call STMSEQ('(la2),2;(la3),3;(la4),20')
cref        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
        call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      case(4)
cold        call STMSEQ('(la2),2;(la3),3;(la4),20')
cref        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
        call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      case(5)
cold        call STMSEQ('(la2),2;(la3),3;(la4),20')
cref         call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
        call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      case(6)
cold        call STMSEQ('(la2),2;(la3),3;(la4),20')
cref        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
        call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
      case(7)
cref        call STMSEQ('(la2),2;(la3),20')
        call STMSEQ('(trp),2;(la3),20')
      case(8)
cref        call STMSEQ('(la2),2;(la3),20')
        call STMSEQ('(trp),2;(la3),20')
      case(9)
cref        call STMSEQ('(la2),2;(la3),20')
        call STMSEQ('(trp),2;(la3),20')
      case(10)
cold        call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref        call STMSEQ('(la2),2;(la3),3;(la4),20')
        call STMSEQ('(trp),2;(la3),3;(la4),20')
      end select
c
      TITLE(1) = 'Chemical Reactor--Bounded Control, Citron'
c
c             define parameters for this case
c
      call CHMRdF(icase)
c
      TITLE(2)(1:5) = 'Case '
      write(title(2)(6:7),'(i2)') icase
      title(2)(8:14) = '   k = '
      write(title(2)(15:17),'(f3.1)') expk
      title(2)(18:25) = '   au = '
      write(title(2)(26:28),'(f3.1)') aupr
      title(2)(29:36) = '   tf = '
      write(title(2)(37:39),'(f3.1)') tfinal
C
      STLBL(0) = 'TIME    Time'
      STLBL(1) = 'XPROD   Chemical product x'
      STLBL(2) = 'YPROD   Chemical product y'
      STLBL(3) = 'ARATE   Rate coefficient a'
c
      NPV = 0
      INIT = 1
      NDE = 2
      NAV = 1
      NCF(1) = NDE
      npf(1) = 0
C
C             Default State and Control Bounds
C
      DO 75 I=-1,1
c
c           ---x-product
c
        YLB(I,1) = -.1d0
        YUB(I,1) = 1.1d0
c
c           ---y-product
c
        YLB(I,2) = -.1d0
        YUB(I,2) = 1.1d0
c
c           ---rate coefficient
c
        YLB(I,3) = alwr
        YUB(I,3) =  aupr
c
 75   CONTINUE
C
C             Initial time and final time (fixed).
C
      Y0(0) = ZERO
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
c
      T1 = tfinal
      Y1(0) = T1
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C             Define initial conditions for state variables, Y.
C             and fix initial values
C
      Y0(1) = xzero
      Y0(2) = yZERO
      DO 20 I=1,nde
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
   20 CONTINUE
C
C             Define terminal conditions for state variables.
C
      call dcopy(nde,y0(1),1,y1(1),1)
C
C             Define initial and final control
C
      Y0(3) = 1.d0
      Y1(3) = 1.d0
C
C             Define final value of y product as objective to maximize.
C
      ITERM(1,1) = 0
      ITERM(2,1) = 1
      ITERM(3,1) = 1
      ITERM(4,1) = 2
      COEF(1) = ONE
C
      MAXMIN = 1
C
C             Load state scale weights for output.
C
      call dfill(nde+nav,one,stskl,1)
C
      RETURN
      END
