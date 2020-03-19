
      function qsubi(i,t)
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
      xsubi = dble(i-1)/dble(nxgrid-1)
      cospix = cos(pi*xsubi)
      pisqr = pi**2
      e2rt = exp(two*rho*t)
c
      term1 = rho*(acoef1 + two*acoef2) + pisqr*(acoef3 + two*acoef4)
      term1 = term1*cospix*exp(rho*t)
c
      term2 = acoef4*pisqr*e2rt
c
      term3 = two*acoef4*pisqr + rho*acoef2
      term3 = term3*e2rt*(cospix)**2
c
      qsubi = term1 - term2 + term3
c
      RETURN
      END
