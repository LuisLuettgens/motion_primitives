

      subroutine kplrde(iphase,eccent,y,ny,p,np,frhs,nfrhs,iferr)
c
c         computes the right hand sides of the algebraic equation
c
c     ******************************************************************
      implicit doubleprecision (a-h,o-z)
c
      parameter  (zero=0.d0,half=0.5d0,one=1.d0)
c
      dimension  y(ny),p(np),frhs(nfrhs)
c
c     ******************************************************************
c
c             compute rhs of algebraic equations
c
      iferr = 0
      anomm = one
c
      frhs(1) = y(1) - eccent*sin(y(1)) - anomm
c
      return
      end
