
      SUBROUTINE MEDIIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,NGRID,
     &                 INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     &                 MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,CUB,
     &                 CLBL,MXTERM,COEF,ITERM,TITLE,IER)
C
C        Minimum Energy Double Integrator Example.  Zhao & Tsiotras, and/or
c        Bryson & Ho,"APPLIED OPTIMAL CONTROL", p 122.
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
      parameter (zero=0.d0,one=1.d0,two=2.d0,sixth=one/6.d0,
     $    quartr=.25d0,half=.5d0)
c
      INCLUDE '../commons/odeprb.cmn'
c
      common /medicm/ smlL,capJ,ngrtot,ipmedi
C
C     ******************************************************************
C
      ipmedi = ITPCLS(IPROB)
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),3;(la4),20')
      call STMSEQ('(trp),2;(la3),3;(la4),20')
c
      smlL = tpseed(iprob)
c
      onem4L = one-4.d0*smlL
      if(smlL.ge.quartr) then
        capJ  = two
      elseif(sixth.le.smlL.and.smlL.lt.quartr) then
        capJ  = two + 6.d0*onem4L**2
      elseif(smlL.lt.sixth) then
        capJ  = 4.d0/(9.d0*smlL)
      endif
C
      TITLE(1) = 'Minimum Energy Double Integrator'
      TITLE(2) = 'smlL =                Jstar = '
      write(title(2)(8:20),'(g10.3)') smlL
      write(title(2)(32:50),'(g18.10)') capJ
C
      NPHS = 1
      ngrtot = ngrid
c
      init(1) = 2
c
      NTERM = 0
      NKON = 0
      BIGBND = 1.D0/HDMCON(5)
C
      STLBL(0) = 'TIME    Time'
      STLBL(1) = 'X       Distance'
      STLBL(2) = 'V       Velocity'
      STLBL(3) = 'U       Control'
C
C     ----NUMBER OF DIFFERENTIAL EQUATIONS (STATES)   
C
      NDE = 2
      NCF(1) = NDE
C
C     ----NUMBER OF ALGEBRAIC VARIABLES (CONTROLS)   
C
      NAV = 1
C
C     ----GUESS FOR INITIAL TIME AND BOUNDARY CONDITION
C
      Y0(0) = 0.D0
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
C
C     ----GUESS FOR FINAL TIME AND BOUNDARY CONDITION
C
      T1 = 1.D0
      Y1(0) = T1
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
C
C     ----DEFINE INITIAL CONDITIONS FOR STATE VARIABLES, Y.
C
C
      Y0(1) = 0.D0
      Y0(2) = 1.D0
C
C     ----FIX THE INITIAL STATES
C
      DO 110 I=1,2
        YLB(-1,I) = Y0(I)
        YUB(-1,I) = Y0(I)
110   CONTINUE
C
C     ----DEFINE FINAL CONDITIONS FOR STATE VARIABLES, Y.
C
      Y1(1) = 0.D0
      Y1(2) = -1.D0
C
C     ----DEFINE INITIAL AND FINAL CONTROL variables.
C
      Y0(3) = -2.D0
      Y1(3) = 2.D0
C
C     ----FIX THE FINAL VALUES FOR STATES
C
      DO 120 I=1,2 
        YLB(1,I) = Y1(I) 
        YUB(1,I) = Y1(I)
120   CONTINUE
c
      if(ipmedi.gt.0) then
c
c     ----interior state bound
c
        yub(0,1) = smlL
c
      else
c
      write(*,*) "HIER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
C Angepasst, damit es kompiliert!
c          CALL PTHCON(NTERM,NKON,NCF,IPHASE,ITERM,MXTERM,COEF,
c     $    CLB,CUB,CLBL,MXPCON,-bigbnd,smlL,one,'XBND',
c     $    'x(t) < L',IERPTH)
c
          if(ierpth.ne.0) then
            print *,'ierpth = ',ierpth
            stop
          endif
c
      endif
C
C     ----DEFINE quadrature AS THE OBJECTIVE TO BE MINIMIZED
C
      MAXMIN = -1
C
      ncf(3) = ncf(3) + 1
      NTERM = NTERM + 1
      ITERM(1,NTERM) = 0
      ITERM(2,NTERM) = iphase
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = -ncf(2)-ncf(3)
      COEF(NTERM) = .5D0
C
      RETURN
      END
