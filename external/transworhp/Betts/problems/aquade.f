

      subroutine aquade(iphase,time,yvec,nyvec,p,np,frhs,nfrhs,iferr)
c
c         computes the right hand sides of the differential equations.
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
      parameter  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
      parameter  (uxmax=2.d0, uzmax=1.d0,csubx=.5d0,csubz=.1d0,
     $            rsubx=.1d0)
c
      dimension  yvec(nyvec),p(np),frhs(nfrhs)
      common /aquacm/ iclass
c
c     ******************************************************************
c
      iferr = 0
c
      arg = -((yvec(1)-csubx)/rsubx)**2 
      if(arg.lt.hdmcon(10)) then
        iferr = 1
        return
      else
        exparg = exp(arg)
      endif
c      
      rx = -uxmax*exparg
      rx = rx*(yvec(1)-csubx)*((yvec(3)-csubz)/csubz)**2
c
      rz = -uzmax*exparg
      rz = rz*((yvec(3)-csubz)/csubz)**2
c
      frhs(1) = cos(yvec(6))*cos(yvec(5))*yvec(7) + rx
      frhs(2) = sin(yvec(6))*cos(yvec(5))*yvec(7)     
      frhs(3) = -sin(yvec(5))*yvec(7) + rz
      frhs(4) = yvec(8) + sin(yvec(4))*tan(yvec(5))*yvec(9)
     $         + cos(yvec(4))*tan(yvec(5))*yvec(10)
      frhs(5) = cos(yvec(4))*yvec(9) - sin(yvec(4))*yvec(10)
      frhs(6) = sin(yvec(4))*yvec(9)/cos(yvec(5))
     $         + cos(yvec(4))*yvec(10)/cos(yvec(5))
c
      call dcopy(4,yvec(11),1,frhs(7),1)
c
c         quadrature objective function
c
      if(iclass.eq.1) then
        frhs(11) = ddot(4,yvec(11),1,yvec(11),1)
      else
        do ii = 11,14
          frhs(ii) = yvec(ii)**2
        enddo
      endif
c
      return
      end
