      
      
      double precision function spyke(t,a,b)
c
      implicit double precision (a-h,o-z)
c
      arg = - b*(t-a)**2
c
      arg = max(arg,hdmcon(10))
c
      spyke = exp(arg)
c
      return
      end      
