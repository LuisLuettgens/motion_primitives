

      subroutine hangde(iphase,t,y,ny,p,np,f,nf,iferr)
c
c         computes the right hand sides of the hang glider
c         differential equations
c
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
c  arguments:
      integer    iphase,ny,np,nf,iferr
      dimension  y(ny),p(np),f(nf)
      dimension  rhs(4)
c
      common /hangcm/ icl
c
      parameter (zero = 0.d0, one=1.d0, two = 2.d0)
      parameter (uamax = 2.5d0, bigr = 100.d0, cd0 = .034d0, 
     $    xk = .069662d0, xmass = 100.d0, sref = 14.d0,
     $    rho = 1.13d0, grav = 9.80665d0)
      parameter (wt = xmass*grav, rhosr2 = rho*sref/two)
c
c     ******************************************************************
c
      if(icl.eq.1) then
        xdist = y(1)
        yalt = y(2)
        vsubx = y(3)
        vsuby = y(4)
        csubl = y(5)
      elseif(icl.eq.2) then
        xdist = y(1)
        yalt = y(2)
        vsubx = y(3)
        vsuby = y(4)
        time = y(5)
        csubl = y(6)
      elseif(icl.eq.3) then
        xdist = t
        yalt = y(1)
        vsubx = y(2)
        vsuby = y(3)
        csubl = y(4)
      endif
c
c       aerodynamic quantities
c
      csubd = cd0 + xk*csubl**2
c
      xrarg = (xdist/bigr - 2.5d0)**2
c
      if(-xrarg.lt.hdmcon(10)) then
        iferr = 1
        return
      else
        usuba = uamax*exp(-xrarg)*(one - xrarg)
      endif
c
      vy = vsuby - usuba
      vrsqr = vsubx**2 + vy**2
c
      vr = sqrt(vrsqr)
c
c       lift and drag
c
      xlift = csubl*rhosr2*vrsqr
      drag = csubd*rhosr2*vrsqr
c
      if(vr.gt.hdmcon(5)) then
        sineta = vy/vr
        coseta = vsubx/vr
      else
        iferr = 2
        return
      endif
c
c             compute state equations.
c
      rhs(1) = vsubx
      rhs(2) = vsuby
      rhs(3) = (-xlift*sineta - drag*coseta)/xmass
      rhs(4) = (xlift*coseta - drag*sineta - wt)/xmass
      if(icl.eq.1) then
        f(1:4) = rhs(1:4)
      elseif(icl.eq.2) then
        f(1:4 ) = time*rhs(1:4)
        f(5) = zero
      elseif(icl.eq.3) then
        f(1) = rhs(2)/vsubx
        f(2) = rhs(3)/vsubx
        f(3) = rhs(4)/vsubx
      endif
c
      iferr = 0
c
      return
      end
