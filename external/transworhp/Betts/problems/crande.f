

      subroutine crande(iphase,time,yvec,nyvec,p,np,frhs,nfrhs,iferr)
c
c         computes the right hand sides of the differential equations.
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
      parameter  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
c
      dimension  yvec(nyvec),p(np),frhs(nfrhs)
      common /crancm/ ccoef
c
c     ******************************************************************
c
      iferr = 0
c
      uone = yvec(7)
      utwo = yvec(8)
c
      frhs(1) = yvec(4)
      frhs(2) = yvec(5)    
      frhs(3) = yvec(6)
      frhs(4) = uone + 17.2656d0*yvec(3)
      frhs(5) = utwo
      if(abs(yvec(2)).lt.hdmcon(5)) then
        iferr = 1
        return
      else
        frhs(6) = -(uone + 27.0756d0*yvec(3) + two*yvec(5)*yvec(6))
     $             /yvec(2)
      endif
c
c         quadrature objective function
c
      frhs(nfrhs) = half*(yvec(3)**2 + yvec(6)**2 + 
     $                 ccoef*(uone**2 + utwo**2))
c
      return
      end
