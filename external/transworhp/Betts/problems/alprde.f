      
      
      
      SUBROUTINE ALPRDE(IPHASE,T,Z,NZ,P,NP,F,NF,IFERR)
      implicit double precision (a-h,o-z)      
 
      
      INTEGER IPHASE,NZ,NP,NF,IFERR,IPATH
      DOUBLE PRECISION T,Z(NZ),P(NP),F(NF),PT312,PT610,PT106
      DOUBLE PRECISION PT154
c
      iferr = 0
      pt312 = spyke(t,3.d0,12.d0)
      pt610 = spyke(t,6.d0,10.d0)
      pt106 = spyke(t,10.d0,6.d0)
      pt154 = spyke(t,15.d0,4.d0)
c
      IPATH = 0
      IPATH = IPATH + 1
      F(IPATH) = -10.d0*Z(1) + Z(5) + Z(6)
      IPATH = IPATH + 1
      F(IPATH) = -2.d0*Z(2) + Z(5) + 2.d0*Z(6)
      IPATH = IPATH + 1
      F(IPATH) = -3.d0*Z(3) + 5.d0*Z(4) + Z(5) - Z(6)
      IPATH = IPATH + 1
      F(IPATH) = 5.d0*Z(3) - 3.d0*Z(4) + Z(5) + 3.d0*Z(6)
      IPATH = IPATH + 1
      F(IPATH) = Z(1)**2 + Z(2)**2 + Z(3)**2 + Z(4)**2 - 3.d0*PT312
     +           -3.d0*PT610 - 3.d0*PT106 - 8.d0*PT154 -.01d0
      IPATH = IPATH + 1
      F(IPATH) = 100.d0*(Z(1)**2 + Z(2)**2 + Z(3)**2 + Z(4)**2)
     +             +  .01d0*(Z(5)**2 + Z(6)**2)
	

      
      RETURN
      END
