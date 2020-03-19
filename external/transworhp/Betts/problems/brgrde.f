

      SUBROUTINE BRGRDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the Burger's equation.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,HALF=0.5D0,one=1.d0)
C
      COMMON /BRGRCM/ EPS,ICL,NDEQ,NGPT
C
      DIMENSION  Y(NY),P(np),F(NF)
C
C     ******************************************************************
C
C             Compute DY.
C
      iferr = 0
      if(icl.eq.1.or.icl.eq.4) then
c
c         standard formulation
c
        call brgreq(iphase,t,y,ny,p,np,f,nf,iferr)
c
      elseif(icl.eq.2) then
c
c         implicit error equidistribution formulation
c
        nym1 = ny - 1
        nfm1 = nf - 1
        time = y(ny)
        call brgreq(iphase,time,y,nym1,pdum,1,f,nfm1,iferr)
c
c         compute arc length monitor function
c
        phi = sqrt(one + ddot(nfm1,f,1,f,1))
c
c         form final differential equation
c
        f(ny) = p(1)/phi
c
c         scale the original set of ode's
c
        call dscal(nym1,f(ny),f,1)
c
      elseif(icl.eq.3) then
c
c         elastic programming mode
c
        F(1) = Y(2) + y(3) - y(4)
        F(2) = Y(1)*Y(2)/EPS + y(5) - y(6)
        f(3) = y(3) + y(4) + y(5) + y(6)
c
      endif
C
      RETURN
      END
