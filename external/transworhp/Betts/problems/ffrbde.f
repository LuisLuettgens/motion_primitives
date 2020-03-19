

      subroutine ffrbde(iphase,time,yvec,nyvec,p,np,frhs,nfrhs,iferr)
c
c         computes the right hand sides of the differential equations.
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
      parameter  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
c
      dimension  yvec(nyvec),p(np),frhs(nfrhs)
      parameter (alfa=.2d0, beta=.2d0)
c
c     ******************************************************************
c
      iferr = 0
c
      uone = yvec(7) - yvec(8)
      utwo = yvec(9) - yvec(10)
      absu1 = yvec(7) + yvec(8)
      absu2 = yvec(9) + yvec(10)
c
      frhs(1) = yvec(4)
      frhs(2) = yvec(5)    
      frhs(3) = yvec(6)
      frhs(4) = (uone + utwo)*cos(yvec(3))
      frhs(5) = (uone + utwo)*sin(yvec(3))
      frhs(6) = alfa*uone - beta*utwo
c
c         quadrature objective function
c
      frhs(nfrhs) = absu1 + absu2
c
      return
      end
