
      subroutine zrmlin(iphase,nphs,method,nstg,ncf,npf,npv,nav,ngrid,
     &                 init,maxmin,mxparm,p0,plb,pub,plbl,
     &                 mxstat,y0,y1,ylb,yub,stskl,stlbl,mxpcon,clb,cub,
     &                 clbl,mxterm,coef,iterm,title,ier)
c
C             Initialize Zermelo's Problem
C
C             REF:  See Example 1, p. 77  in Bryson & Ho
c                   "applied optimal control"
c
c     ******************************************************************
      implicit doubleprecision (a-h,o-z)
c
c  arguments:
      integer    iphase,nphs,method,nstg,ncf(3),npf(2),npv,nav,ngrid,
     &           init,maxmin,mxparm,mxstat,mxpcon,mxterm,
     &           iterm(4,mxterm),ier
      dimension  p0(mxparm),plb(mxparm),pub(mxparm),y0(0:mxstat),
     &           y1(0:mxstat),ylb(-1:1,0:mxstat),yub(-1:1,0:mxstat),
     &           stskl(0:mxstat),clb(mxpcon),cub(mxpcon),coef(mxterm)
      character  title(3)*60,plbl(mxparm)*80,stlbl(0:mxstat)*80,
     &           clbl(0:mxpcon)*80
c
      INCLUDE '../commons/odeprb.cmn'
c
c     ******************************************************************
c
C             Define Labels for Problem
C
      TITLE(1) = "Zermelo's Problem"
      TITLE(2) = 'Linear Guess'
      STLBL(0) = 'TIME    Time'
      STLBL(1) = 'XPOS    X-Position'
      STLBL(2) = 'YPOS    Y-Position'
      STLBL(3) = 'THETA   Heading Angle, Theta'
      PLBL(1)  = 'TFINAL  Final time'
c
      nterm = 0
      nkon = 0
c
c             number of differential equations, and algebraic variables
c
      nde = 2
      nav = 1
      ncf(1) = nde
c
c             initial time (fixed) and guess for final time.
c
      y0(0) = 0.d0
      ylb(-1,0) = y0(0)
      yub(-1,0) = y0(0)
      t1 = 5.d0
      y1(0) = t1
      ylb(1,0) = 0.d0
c
c             define initial and final conditions for state variables, y.
c             also fix the initial and final values.
c
      y0(1) = 3.5d0
      y0(2) = -1.8d0
      y1(1) = 0.d0
      y1(2) = 0.d0
      do 20 i=1,nde
        ylb(-1,i) = y0(i)
        yub(-1,i) = y0(i)
        ylb(1,i) = y1(i)
        yub(1,i) = y1(i)
   20 continue
c
c             control guess (2*pi/3)
c
      y0(3) = 3.14d0*2.d0/3.d0
      y1(3) = 3.14d0*2.d0/3.d0
c
c             display the output heading angle in degrees
c
      cr2d =  57.29577951308232d0
      stskl(3) = cr2d
c
c             define final time as objective.
c
      nterm = nterm + 1
      iterm(1,nterm) = 0
      iterm(2,nterm) = 1
      iterm(3,nterm) = 1
      iterm(4,nterm) = 0
      coef(nterm) = 1.d0
c     
      maxmin = -1
c
      return
      end
