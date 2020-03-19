
      subroutine zrmlde(iphase,t,y,ny,p,np,f,nf,iferr)
c
c         computes the right hand sides of the 
c         differential equations for Zermelo's problem.
c
c     ******************************************************************
      implicit doubleprecision (a-h,o-z)
c
c  arguments:
      integer    iphase,ny,np,nf,iferr
      dimension  y(ny),p(np),f(nf)
c
c     ******************************************************************
c
      iferr = 0
c
c             compute right hand sides
c
      theta = y(3)
c
      vel = 1.d0
      constk = -1.d0
      xdot = vel*cos(theta) + constk*y(2)
      ydot = vel*sin(theta)
c
      f(1) = xdot
      f(2) = ydot
c
      return
      end
