

      SUBROUTINE TB2SDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the two-strain tuberculosis
c         model
C
C     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ARGUMENTS:
C
      INTEGER    IPHASE,NY,NP,NF,IFERR
      DIMENSION  Y(NY),P(NP),F(NF)
      parameter  (beta1 = 13.d0, beta2 = 13.d0)
      parameter  (xmu = .0143d0, done  = 0.d0 , dtwo   = 0.d0  )
      parameter  (xk1 = .5d0,    xk2   = 1.d0 , rone   = 2.d0  )
      parameter  (rtwo = 1.d0,   smlp  = .4d0 , smlq   = .1d0  )
      parameter  (capN = 30.d3,  bstar  = .029d0 )
      parameter  (Bone  = 50.d0, Btwo   = 500.d0)
      parameter  (capLam = xmu*capN)
      parameter  (Szero  = 76.d0*capN/120.d0)
      parameter  (capL10 = 36.d0*capN/120.d0)
      parameter  (capI10 =  4.d0*capN/120.d0)
      parameter  (capL20 =  2.d0*capN/120.d0)
      parameter  (capI20 =  1.d0*capN/120.d0)
      parameter  (capT0  =  1.d0*capN/120.d0)
      parameter  (zero = 0.d0, one = 1.d0)
C
C     ******************************************************************
C
C     ----SET FUNCTION ERROR FLAG.
C
      IFERR = 0
C
C     ----Unload dynamic variables into local quantities
C
      capS = y(1)
      capL1 = y(2)
      capI1 = y(3)
      capL2 = y(4)
      capI2 = y(5)
      capT  = y(6)
c
      uone  = y(7)
      utwo  = y(8)
C
C     ----COMPUTE STATE EQUATIONS.
C
      dotS  = capLam - beta1*capS*capI1/capN 
     $      - bstar*capS*capI2/capN - xmu*capS
c
      dotL1 = beta1*capS*capI1/capN - (xmu + xk1)*capL1 
     $      - uone*rone*capL1 + (one - utwo)*smlp*rtwo*capI1
     $      + beta2*capT*capI1/capN - bstar*capL1*capI2/capN
c
      dotI1 = xk1*capL1 - (xmu + done)*capI1 - rtwo*capI1
c
      dotL2 = (one - utwo)*smlq*rtwo*capI1 - (xmu + xk2)*capL2
     $      + bstar*(capS + capL1 + capT)*capI2/capN
c
      dotI2 = xk2*capL2 - (xmu + dtwo)*capI2  
c
      dotT  = uone*rone*capL1 - xmu*capT 
     $      + (one - (one - utwo)*(smlp + smlq))*rtwo*capI1 
     $      - beta2*capT*capI1/capN - bstar*capT*capI2/capN            
C
      F(1) = dotS
      F(2) = dotL1
      F(3) = dotI1           
      F(4) = dotL2         
      F(5) = dotI2           
      F(6) = dotT     
c
c     ----quadrature integrand
c
      f(7) = capL2 + capI2 + Bone*uone**2/2.d0 + Btwo*utwo**2/2.d0
C
      RETURN
      END
