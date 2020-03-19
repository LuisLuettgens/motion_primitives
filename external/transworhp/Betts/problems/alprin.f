      
      SUBROUTINE ALPRIN(IPHASE,NPHS,METHOD,NSTG,NCF,NPF,NPV,NAV,
     +                  NGRID,INIT,MAXMIN,MXPARM,P0,PLB,PUB,PLBL,
     +                  MXSTAT,Y0,Y1,YLB,YUB,STSKL,STLBL,MXPCON,CLB,
     +                  CUB,CLBL,MXTERM,COEF,ITERM,TITLE,IER)
     
     
      implicit double precision (a-h,o-z)      
          
      INTEGER IPHASE,NPHS,METHOD,NSTG,NCF(3),NPF(2),NPV,NAV,NGRID,
     +        INIT,MAXMIN,MXPARM,MXSTAT,MXPCON,MXTERM,ITERM(4,MXTERM),
     +        IER,NTERM,NPATH
     
      DOUBLE PRECISION P0(MXPARM),PLB(MXPARM),PUB(MXPARM),Y0(0:MXSTAT),
     +          Y1(0:MXSTAT),YLB(-1:1,0:MXSTAT),YUB(-1:1,0:MXSTAT),
     +          STSKL(0:MXSTAT+MXPARM,2),CLB(MXPCON),CUB(MXPCON),
     +          COEF(MXTERM)
     
      CHARACTER PLBL(MXPARM+2)*80,STLBL(0:MXSTAT)*80,
     +          CLBL(0:MXPCON)*80,TITLE(3)*60
     
      INCLUDE '../commons/odeprb.cmn'
c
      parameter (zero=0.d0, one=1.0d0)
c
      icls = itpcls(iprob)
      method = itpmet(iprob)
      nstg = itpstg(iprob)
      ngrid = itpgrd(iprob)
c
cold      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
cref      call STMSEQ('(la2),2;(la3),20')
      call STMSEQ('(trp),2;(la3),20')
c
      NCF(1) = 4
      NAV = 2
c      
      TITLE(1) = 'Alp Rider (Basic)'
      TITLE(2) = 'Minimum Integral with Path Inequality'
c
      INIT = 1
      BIGBND = ONE/HDMCON(5)
   
      
      Y0(0) = 0.d0
      Y0(1) = 2.d0
      Y0(2) = 1.d0
      Y0(3) = 2.d0
      Y0(4) = 1.d0
      Y0(5) = 0.d0
      Y0(6) = 0.d0
      
      Y1(0) = 20.d0
      Y1(1) = 2.d0
      Y1(2) = 3.d0
      Y1(3) = 1.d0
      Y1(4) = -2.d0
      Y1(5) = 0.d0
      Y1(6) = 0.d0
   
      
      
      YLB(-1,0) = Y0(0)
      YUB(-1,0) = Y0(0)
      YLB(1,0) = Y1(0)
      YUB(1,0) = Y1(0)
      YLB(-1,1) = Y0(1)
      YUB(-1,1) = Y0(1)
      YLB(-1,2) = Y0(2)
      YUB(-1,2) = Y0(2)
      YLB(1,1) = Y1(1)
      YUB(1,1) = Y1(1)
      YLB(1,2) = Y1(2)
      YUB(1,2) = Y1(2)
      YLB(-1,3) = Y0(3)
      YUB(-1,3) = Y0(3)
      YLB(1,3) = Y1(3)
      YUB(1,3) = Y1(3)
      YLB(-1,4) = Y0(4)
      YUB(-1,4) = Y0(4)
      YLB(1,4) = Y1(4)
      YUB(1,4) = Y1(4)
     
      NTERM = 0
      NPATH = 0
      NCF(2) = 0
      NTERM = NTERM + 1
      NPATH = NPATH+1
      NCF(2) = NCF(2) + 1
      ITERM(1,NTERM) = 1
      ITERM(2,NTERM) = 1
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = -NPATH
      COEF(NTERM) = 1.d0
      CLB(1) = 0.d0
      clBL(1) = 'PATHCN  Terrain Path Constraint'
     

       
      NCF(3) = 1
      NTERM = NTERM + 1
      NPATH = NPATH+1
      ITERM(1,NTERM) = 0
      ITERM(2,NTERM) = 1
      ITERM(3,NTERM) = 0
      ITERM(4,NTERM) = -NPATH
      COEF(NTERM) = 1.d0
      
      RETURN
      END
