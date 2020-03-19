

      SUBROUTINE TRANDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         COMPUTES THE RIGHT HAND SIDES OF THE train
c         DIFFERENTIAL EQUATIONS.
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ARGUMENTS:
C
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
C
      dimension zhill(2),shill(3)
c
      PARAMETER (acoef = .3d0, bcoef = .14d0, ccoef = .16d0,
     $     eps = .05d0)
      data (zhill(ii),ii=1,2) / 2.d0, 4.d0 /
      data (shill(ii),ii=1,3) / 2.d0, 0.d0, -2.d0 /
C
C     ******************************************************************
C
C     ----SET FUNCTION ERROR FLAG.
C
      IFERR = 0
c
      x = y(1)
      v = y(2)
      ua = y(3)
      ub = y(4)
C
C     ----COMPUTE hill function.
C
      pi = hdmcon(12)
      hx = 0.d0
      do j = 1,2
        term = (shill(j+1)-shill(j))*atan2(x-zhill(j),eps)/pi
        hx = hx + term
      enddo
C
C     ----COMPUTE STATE EQUATIONS.
C
      F(1) = v
      F(2) = hx - (acoef + bcoef*v + ccoef*v**2) + ua - ub
c
c     ----define quadrature function
c
      tol = 1.d-3
      F(3) = ua*v + tol*(ua**2 + ub**2)
C     
C
      RETURN
      END
