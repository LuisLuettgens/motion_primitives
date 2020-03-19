
      subroutine hangin(iphase,nphs,method,nstg,ncf,npf,npv,nav,ngrid,
     &                 init,maxmin,mxparm,p0,plb,pub,plbl,
     &                 mxstat,y0,y1,ylb,yub,stskl,stlbl,mxpcon,clb,cub,
     &                 clbl,mxterm,coef,iterm,title,ier)
c
c             initialize the hang glider problem.
c             ref.  "combining direct and indirect methods in optimal
c             control:  range maximization of a hang glider,"  r. bulirsch,
c             e. nerz, h.j. pesch, o. von stryk,  rept. no. 313, 1991,
c             mathematisches institut, technische universitat munchen, 
c             arcisstrasse 21, d-8000, muenchen 2
c
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
      parameter (zero=0.d0,one=1.d0,two=2.d0)
c
c  arguments:
      integer    iphase,nphs,method,nstg,ncf(3),npf(2),npv,nav,ngrid,
     &           init(2),maxmin,mxparm,mxstat,mxpcon,mxterm,
     &           iterm(4,mxterm),ier
      dimension  p0(mxparm),plb(mxparm),pub(mxparm),y0(0:mxstat),
     &           y1(0:mxstat),ylb(-1:1,0:mxstat),yub(-1:1,0:mxstat),
     &           stskl(0:mxstat+mxparm,2),clb(mxpcon),cub(mxpcon),
     &           coef(mxterm)
      character  title(3)*60,plbl(mxparm)*80,stlbl(0:mxstat)*80,
     &           clbl(0:mxpcon)*80
c
      include '../commons/odeprb.cmn'
c
      common /hangcm/ icl
c
c     ******************************************************************
c
c             obtain number of grid points for problem.
c
      nphs = 1
      icl = itpcls(iprob)
      ngrid = itpgrd(iprob)
      method = itpmet(iprob)
      nstg = itpstg(iprob)
c
      call STMSEQ('(trp),1;(la4),20')
c
      title(1) = 'maximum range flight of a hang glider'
c
      npv = 0
      init(1) = 1
      if (icl.eq.1) then
          title(2) = ' '
          nde = 4
          nav = 1
      elseif (icl.eq.2) then
          title(2) = 'Formulation in Standard Form'
          nde = 5
          nav = 1
      elseif (icl.eq.3) then
          title(2) = 'Independent Range Formulation'
          nde = 3
          nav = 1
      else
          ier = -11
          return
      endif
      ncf(1) = nde
c
      if(icl.eq.1) then
c
c             initial time (fixed) and guess for final time.
c
        y0(0) = zero
        ylb(-1,0) = y0(0)
        yub(-1,0) = y0(0)
        t1 = 100.d0
        y1(0) = t1
        ylb(1,0) = zero
        yub(1,0) = 1.1d0*t1
c
c             define initial conditions for state variables, y.
c                 y(1) = horizontal distance (x in meters)
c                 y(2) = altitude (y in meters)
c                 y(3) = horizontal velocity. (meters/sec)
c                 y(4) = vertical velocity. (meters/sec)
c
c             declare initial values as fixed.
c
        y0(1) = zero
        y0(2) = 1000.d0
        y0(3) = 13.227567500d0
        y0(4) = -1.2875005200d0
        do 20 i=1,4
          ylb(-1,i) = y0(i)
          yub(-1,i) = y0(i)
   20   continue
c
c             define terminal conditions for state variables.
c
        y1(1) = 1250.d0
        y1(2) = 900.d0
        y1(3) = y0(3)
        y1(4) = y0(4)
c
c
c               define initial and final control (csubl)
c
        y0(5) = one
        y1(5) = y0(5)
c
c             fix final states 2-4.
c
        do 70 i=2,4
          ylb(1,i) = y1(i)
          yub(1,i) = y1(i)
   70   continue
c
c             bound control magnitude on phase by 90.
c
        do 75 i=-1,1
          ylb(i,5) = zero
          yub(i,5) = 1.4d0
   75   continue
c
c             define final range as objective (maximize)
c
        iterm(1,1) = 0
        iterm(2,1) = 1
        iterm(3,1) = 1
        iterm(4,1) = 1
        coef(1) = one
c
        maxmin = 1
c
      elseif(icl.eq.2) then
c
c             initial and final tau fixed and guess for final time.
c
        y0(0) = zero
        ylb(-1,0) = y0(0)
        yub(-1,0) = y0(0)
        t1 = 100.d0
        y1(0) = one
        ylb(1,0) = one
        yub(1,0) = one
c
c             define initial conditions for state variables, y.
c                 y(1) = horizontal distance (x in meters)
c                 y(2) = altitude (y in meters)
c                 y(3) = horizontal velocity. (meters/sec)
c                 y(4) = vertical velocity. (meters/sec)
c                 y(5) = time (sec)
c
c             declare initial values as fixed.
c
        y0(1) = zero
        y0(2) = 1000.d0
        y0(3) = 13.227567500d0
        y0(4) = -1.2875005200d0
        y0(5) = t1
        do 120 i=1,4
          ylb(-1,i) = y0(i)
          yub(-1,i) = y0(i)
  120   continue
c
c             define terminal conditions for state variables.
c
        y1(1) = 1250.d0
        y1(2) = 900.d0
        y1(3) = y0(3)
        y1(4) = y0(4)
        y1(5) = t1
c
c
c               define initial and final control (csubl)
c
        y0(6) = one
        y1(6) = y0(6)
c
c             fix final states 2-4.
c
        do 170 i=2,4
          ylb(1,i) = y1(i)
          yub(1,i) = y1(i)
  170   continue
c
c             bound on final time
c
        ylb(1,5) = zero
        yub(1,5) = 1.1d0*t1
c
c             bound control magnitude on phase
c
        do 175 i=-1,1
          ylb(i,6) = zero
          yub(i,6) = 1.4d0
  175   continue
c
c             define final range as objective (maximize)
c
        iterm(1,1) = 0
        iterm(2,1) = 1
        iterm(3,1) = 1
        iterm(4,1) = 1
        coef(1) = one
c
        maxmin = 1
c
      elseif(icl.eq.3) then
c
c             initial range (fixed) and guess for final range
c
        y0(0) = zero
        ylb(-1,0) = y0(0)
        yub(-1,0) = y0(0)
        y1(0) = 1250.d0
        ylb(1,0) = zero
        yub(1,0) = 1500.d0
c
c             define initial conditions for state variables, y.
c                 y(1) = altitude (y in meters)
c                 y(2) = horizontal velocity. (meters/sec)
c                 y(3) = vertical velocity. (meters/sec)
c
c             declare initial values as fixed.
c
        y0(1) = 1000.d0
        y0(2) = 13.227567500d0
        y0(3) = -1.2875005200d0
c
        ylb(-1,1:3) = y0(1:3)
        yub(-1,1:3) = y0(1:3)
c
c             define terminal conditions for state variables.
c
        y1(1) = 900.d0
        y1(2) = y0(2)
        y1(3) = y0(3)
c
c
c               define initial and final control (csubl)
c
        y0(4) = one
        y1(4) = y0(4)
c
c             fix final states 1-3.
c
        ylb(1,1:3) = y1(1:3)
        yub(1,1:3) = y1(1:3)
c
c             bound control magnitude on phase by 90.
c
        ylb(-1:1,4) = zero
        yub(-1:1,4) = 1.4d0
c
c             define final range as objective (maximize)
c
        iterm(1,1) = 0
        iterm(2,1) = 1
        iterm(3,1) = 1
        iterm(4,1) = 0
        coef(1) = one
c
        maxmin = 1
c
      endif
c
      return
      end
