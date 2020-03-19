

      subroutine lnhtde(iphase,time,yvec,nyvec,p,np,frhs,nfrhs,iferr)
c
c         computes the right hand sides of the differential equations.
c
c     ******************************************************************
      implicit double precision (a-h,o-z)
c
      parameter  (zero=0.d0,half=0.5d0,one=1.d0,two=2.d0)
c
      dimension  yvec(nyvec),p(np),frhs(nfrhs)
c
      common /lnhtcm/ icls
c
c     ******************************************************************
c
      iferr = 0
c
      zmte = 2.0d-2
      zmbe = 2.4d-1
      zmvau = 2.4d0
      zk1 = 2.4d-5
      zk2 = 3.0d-3
      rcon = 3.0d-2
      tmax = 1500.0d0
      scon = 10.0d0
      zN = 1200.0d0
c
      ucntrl = yvec(nyvec) 
c
      FRHS(1) = scon/(1.0d0+yvec(4)) - zmte*yvec(1) 
     $   + rcon*yvec(1)*(1.0d0 - (yvec(1) + yvec(2) + yvec(3))/tmax) 
     $   - zk1*yvec(4)*yvec(1) 
      FRHS(2) = zk1*yvec(4)*yvec(1) - zmte*yvec(2) - zk2*yvec(2)
      FRHS(3) = zk2*yvec(2) - zmbe*yvec(3)
      FRHS(4) = zN*ucntrl*zmbe*yvec(3) - zk1*yvec(4)*yvec(1) 
     $   - zmvau*yvec(4)
      FRHS(5) = - yvec(1) + 50.0d0*(1.0d0 - ucntrl)**2
      if(icls.eq.1) frhs(5) = 1.d-5*frhs(5)
c
      return
      end
