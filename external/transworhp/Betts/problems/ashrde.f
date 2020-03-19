

      SUBROUTINE ASHRDE(IPHASE,Xvar,Y,NY,P,NP,Frhs,NFrhs,IFERR)
C
C         Computes the right hand sides of the differential equations.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,HALF=0.5D0,one=1.d0,two=2.d0)
C
      COMMON /ASHRCM/ EPS,ICL,NDEQ,NGPT
C
      DIMENSION  Y(NY),P(np),Frhs(NFrhs)
C
C     ******************************************************************
C
C             Compute DY.
C
      iferr = 0
      pi = hdmcon(12)
      pix = pi*xvar
c
      if(icl.lt.3) then
c
        rhs = xvar*y(2) + eps*(pi**2)*cos(pix) + pix*sin(pix)
c
        frhs(1) = y(2)
        frhs(2) = -rhs/eps
c
      elseif(icl.eq.3) then
c
        frhs(1) = y(2)/sqrt(pi*eps)
        frhs(2) = -two*xvar*y(2)/eps
c
      elseif(icl.eq.4) then
c
        frhs(1) = y(2)
        frhs(2) = y(1) - 999.999d0*y(2)
c
      elseif(icl.eq.5) then
c
        frhs(1) = 1.d0 + y(2)*y(1)**2 - 4.d0*y(1)
        frhs(2) = 3.d0*y(1) - y(2)*y(1)**2
c
      elseif(icl.eq.6) then
c
        frhs(1) = 1.d0 + y(2)*y(1)**2 - 4.d0*y(1) 
        frhs(2) = 3.d0*y(1) - y(2)*y(1)**2
c
        frhs(1) = frhs(1) - y(3) + y(4)
        frhs(2) = frhs(2) - y(5) + y(6)
        frhs(3) = y(3) + y(4) + y(5) + y(6)
c
      endif
C
      return
      end
