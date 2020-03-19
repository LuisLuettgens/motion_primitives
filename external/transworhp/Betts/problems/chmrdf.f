


      SUBROUTINE CHMRdF(icase)
C
C         Define the constants for case icase 
C
C     ******************************************************************
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
C
      PARAMETER  (ZERO=0.D0,ONE=1.D0,two=2.d0)
C
C  Arguments:
C
C  Local:
      COMMON /CHMRCM/ alwr,aupr,tfinal,xzero,yzero,rhocof,expk
C
C     ******************************************************************
C
      alwr = .1d0
      xzero = 1.d0
      yzero = .01d0
      rhocof = 2.5d0
      if(icase.eq.1) then
        aupr = .5d0
        tfinal = 2.d0
        expk = 1.5d0
      elseif(icase.eq.2) then
        aupr = .5d0
        tfinal = 4.d0
        expk = 1.5d0
      elseif(icase.eq.3) then
        aupr = .5d0
        tfinal = 8.d0
        expk = 1.5d0
      elseif(icase.eq.4) then
        aupr = .2d0
        tfinal = 2.d0
        expk = 1.5d0
      elseif(icase.eq.5) then
        aupr = .3d0
        tfinal = 2.d0
        expk = 1.5d0
      elseif(icase.eq.6) then
        aupr = .4d0
        tfinal = 2.d0
        expk = 1.5d0
      elseif(icase.eq.7) then
        alwr = .01d0
        aupr = 8.d0
        tfinal = 2.d0
        expk = 1.5d0
      elseif(icase.eq.8) then
        alwr = .01d0
        aupr = 8.d0
        tfinal = 4.d0
        expk = 1.5d0
      elseif(icase.eq.9) then
        alwr = .01d0
        aupr = 8.d0
        tfinal = 8.d0
        expk = 1.5d0
      elseif(icase.eq.10) then
        aupr = .5d0
        tfinal = 2.d0
        expk = .5d0
      else
        print *,'ivalid icase =',icase
        stop
      endif
C
      RETURN
      END
