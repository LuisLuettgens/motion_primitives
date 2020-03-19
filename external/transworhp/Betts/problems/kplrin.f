

      subroutine kplrin(iphase,nphs,method,nstg,ncf,npf,npv,nav,ngrid,
     &                 init,maxmin,mxparm,p0,plb,pub,plbl,
     &                 mxstat,y0,y1,ylb,yub,stskl,stlbl,mxpcon,clb,cub,
     &                 clbl,mxterm,coef,iterm,title,ier)
c
c             initialize kepler's equation problem
c
c
c     ******************************************************************
      implicit doubleprecision (a-h,o-z)
c
      parameter (zero=0.d0,one=1.d0,oneep1=1.d+01,oneep3=1.d+03)
c
c  arguments:
      integer    iphase,nphs,method,nstg,ncf(3),npf(2),npv,nav,ngrid,
     &           init,maxmin,mxparm,mxstat,mxpcon,mxterm,
     &           iterm(4,mxterm),ier
      dimension  p0(mxparm),plb(mxparm),pub(mxparm),y0(0:mxstat),
     &           y1(0:mxstat),ylb(-1:1,0:mxstat),yub(-1:1,0:mxstat),
     &           stskl(0:mxstat+mxparm,2),clb(mxpcon),cub(mxpcon),
     &           coef(mxterm)
      character  title(3)*60,plbl(mxparm)*80,stlbl(0:mxstat)*80,
     &           clbl(0:mxpcon)*80
c
c  local:
c
      include '../commons/odeprb.cmn'
c
c     ******************************************************************
c
c             obtain number of grid points for problem.
c
      nphs = 1
      NGRID = ITPGRD(IPROB)
      METHOD = ITPMET(IPROB)
      NSTG = ITPSTG(IPROB)
c
cold      call STMSEQ('(la2),2;(la3),20')
cold      call STMSEQ('(la2),2;(la3),3;(la4),20')
cref      call STMSEQ('(la2),-1;(la3),-2;(la4),-3;(la5),-20')
      call STMSEQ('(trp),-1;(la3),-2;(la4),-3;(la5),-20')
c
      TITLE(1) = 'Kepler''s equation'
c
      TITLE(2) = 'Algebraic Equation only (no diff. eq.)'
      npv = 0
      nav = 1
      nde = 0
      ncf(1) = 0
      ncf(2) = 1
c
c             initial eccentricity and final eccentricity fixed.
c
      y0(0) = zero
      ylb(-1,0) = y0(0)
      yub(-1,0) = y0(0)
      y1(0) = .9d0
      ylb(1,0) = y1(0)
      yub(1,0) = y1(0)
      STLBL(0) = 'ECCENT  Eccentricity'
c
c             define terminal conditions for state variables (x).
c
      y0(1) = zero
      y1(1) = zero
      STLBL(1) = 'ECCAN   Eccentric Anomoly'
c
c         default constraints 
c
      nterm =0
      nkon = 0
c
c             define constant value path constraint
c
      nterm = nterm + 1
      nkon = nkon + 1
      clbl(nkon) = 'KEPLER  Keplers equation = 0'
      clb(nkon) = zero
      cub(nkon) = zero
      iterm(1,nterm) = nkon
      iterm(2,nterm) = iphase
      iterm(3,nterm) = 0
      iterm(4,nterm) = -1
      coef(nterm) = one
c
      maxmin = 0
c
      init = 1
c
c             load state scale weights for output.
c
      do 150 i=1,1+nde+nav
        stskl(i,1) = one
  150 continue
c
      ier = 0
c
      return
      end
