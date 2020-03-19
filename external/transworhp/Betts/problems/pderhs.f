
      function pderhs(i,t,yip1,yi,yim1)
c
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /HEATC1/ TFINAL,VUPR,CAPH,VLWR,QSUBA,GAMMA,HSTEP,
     &             PI,ACOEF1,ACOEF2,ACOEF3,ACOEF4,GCOEF,RHO,GAMRLX
      COMMON /HEATC2/ ICL,NXGRID
C
      PARAMETER (ZERO = 0.D0, POINT5 = .5D0, ONE = 1.D0, two = 2.d0)
C
      a1pa2y = acoef1 + acoef2*yi
      a3pa4y = acoef3 + acoef4*yi
      ysubxx = (yip1 - two*yi + yim1)/(hstep**2)
      ysubx  = (yip1 - yim1)/(two*hstep)
c
      qit = qsubi(i,t)
c
      pderhs = (qit + a3pa4y*ysubxx + acoef4*ysubx**2)/a1pa2y
c
      RETURN
      END
