
      SUBROUTINE BRGrEQ(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         Computes the right hand sides of the Burger's equation.
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,HALF=0.5D0)
C
      COMMON /BRGRCM/ EPS,ICL,NDEQ,NGPT
C
      DIMENSION  Y(NY),P(np),F(NF)
C
C     ******************************************************************
C
C             Compute DY.
C
      F(1) = Y(2)
      F(2) = Y(1)*Y(2)/EPS
      IFERR = 0
C
      RETURN
      END
