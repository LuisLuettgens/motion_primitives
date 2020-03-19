

      subroutine jshide(iphase,time,yvec,nyvec,p,np,frhs,nfrhs,iferr)
c
c         computes the right hand sides of the differential equations.
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
      parameter  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
c
      dimension  yvec(nyvec),p(np),frhs(nfrhs)
      parameter (s1 = 2.d0, s2 = 1.5d0, xmu = .002d0, conk = 2.5d-4,
     $   ccon = .007d0, gcon = 30.d0, b1 = 14.d0, b2 = 1.d0, a1 = 2.5d5,
     $   a2 = 75.d0)
c
      common /jshicm/ icls
c
c     ******************************************************************
c
      iferr = 0
c
      t = yvec(1)
      v = yvec(2)
      if(icls.eq.1) then
        u1 = yvec(3)
        u2 = yvec(4)
      elseif(icls.eq.2) then
        quad = yvec(3)
        u1 = yvec(4)
        u2 = yvec(5)
      endif
c
      frhs(1) = s1 - (s2*v)/(b1 + v) - xmu*t - conk*v*t + u1*t
      frhs(2) = (gcon*(one-u2)*v)/(b2 + v) - ccon*v*t   
c
c         quadrature objective function
c
      frhs(nfrhs) = t - (a1*u1**2 + a2*u2**2)
c
      return
      end
