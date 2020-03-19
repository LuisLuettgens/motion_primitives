

      SUBROUTINE BRACDE(IPHASE,T,Y,NY,P,NP,F,NF,IFERR)
C
C         COMPUTES DE FOR BRACHISTOCHRONE
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,HALF=0.5D0)
      PARAMETER (TANTHT=.5D0)
C
      DIMENSION  Y(NY),F(NF)
      COMMON /BRACCM/ HBRAC,ICLBRC,NGRDBR
C
C     ******************************************************************
C
C             COMPUTE DY.
C
      IFERR = 0
      GRAV  = 32.174D0
C
      F(1) = Y(3)*COS(Y(4))
      F(2) = Y(3)*SIN(Y(4))
      F(3) = GRAV*SIN(Y(4))
C
      IF(ICLBRC.EQ.2) THEN
        F(4) = Y(2) - Y(1)*TANTHT - HBRAC
      ELSEIF(ICLBRC.EQ.3) THEN
        F(4) = F(2) - F(1)*TANTHT 
      ENDIF
C
      RETURN
      END
