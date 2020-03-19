

      SUBROUTINE HEATDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the heat dynamics system
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Arguments:
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
C
C
      COMMON /HEATC1/ TFINAL,VUPR,CAPH,VLWR,QSUBA,GAMMA,HSTEP,
     &             PI,ACOEF1,ACOEF2,ACOEF3,ACOEF4,GCOEF,RHO,GAMRLX
      COMMON /HEATC2/ ICL,NXGRID
C
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0)
C
C     ******************************************************************
c
      IFERR = 0
C
      if(icl.eq.1) then
C
        LCNTRL = NXGRID + 1
        ONEHSQ = ONE/(HSTEP**2)
C
C             TEMPERATURE PDE'S
C
        F(1) = ONEHSQ*(Y(2) - TWO*Y(1) + Y(LCNTRL+2))
        DO 110 I = 2,NXGRID-1
          F(I) = ONEHSQ*(Y(I+1) - TWO*Y(I) + Y(I-1))
 110    CONTINUE
        F(NXGRID) = ONEHSQ*(Y(LCNTRL+3) - TWO*Y(NXGRID) + Y(NXGRID-1))
C
C             BOUNDARY CONDITION ODE
C
        F(NXGRID+1) = (Y(LCNTRL+1) - Y(NXGRID+1))/GAMMA
C
C             BOUNDARY CONDTION AT X = 0
C
        F(NXGRID+2) = CAPH*(Y(1) - Y(NXGRID+1)) 
     $              - (Y(2) - Y(LCNTRL+2))/(two*HSTEP)
C
C             BOUNDARY CONDTION AT X = 1
C
        F(NXGRID+3) = (Y(LCNTRL+3) - Y(NXGRID-1))/(two*HSTEP)
C
      else
C
        LCNTRL = NXGRID 
C
C             TEMPERATURE PDE'S
C
        F(1) = pderhs(1,t,y(2),y(1),y(lcntrl+2))
        DO 120 I = 2,NXGRID-1
          F(I) = pderhs(i,t,Y(I+1),Y(I),Y(I-1))
 120    CONTINUE
        F(NXGRID) = pderhs(nxgrid,t,Y(lcntrl+3),Y(Nxgrid),Y(Nxgrid-1))
C
C             BOUNDARY CONDTION AT X = 0
C
        F(NXGRID+1) = gcoef*(Y(1) - Y(lcntrl+1)) 
     $       - (acoef3 + acoef4*y(1))*(Y(2) - Y(LCNTRL+2))/(two*HSTEP)
C
C             BOUNDARY CONDTION AT X = 1
C
        F(NXGRID+2) = (acoef3 + acoef4*y(nxgrid))
     $                *(Y(LCNTRL+3) - Y(NXGRID-1))/(two*HSTEP)
C
C             temperature deviation quadrature integrand
C
        ysubd = two - exp(rho*t)
        F(NXGRID+3) = ((y(nxgrid) - ysubd)**2 
     $                + gamrlx*y(lcntrl+1)**2)/two
C
      endif
C
      RETURN
      END
