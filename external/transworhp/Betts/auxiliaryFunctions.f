      subroutine STMSEQ(seq)
C
      character seq(1)*100
C
      return
      end

C
C    *****************************************
C
C     Konstanten aus PDF
      double precision function hdmcon(k)
C
      integer k
C
      select case (k)
      case (1)
        hdmcon = 1.7976000000000000d+308
C        hdmcon = 1.797693134862316D+308
      case (2)
        hdmcon = 4.494232837155790d+307
      case (3)
        hdmcon = 1.7976000000000000d+308
C        hdmcon = 1.797693134862316d+308
      case (4)
        hdmcon = 2.225073858507201d-308
      case (5)
        hdmcon = 2.220446049250313d-016
      case (6)
        hdmcon = 2.220446049250313d-016
      case (7)
        hdmcon = 2.00000000000000
      case (8)
        hdmcon = 53.0000000000000
      case (9)
        hdmcon = 709.000000000000
      case (10)
        hdmcon = -708.000000000000
      case (11)
        hdmcon = 2147483647.00000
      case (12)
        hdmcon = 3.14159265358979
      case (13)
        hdmcon = 2.71828182845905
      case (14)
        hdmcon = 0.577215664901533
      case (15)
        hdmcon = 1.745329251994330d-002
      case (16)
        hdmcon = 57.2957795130823
      case default
        write(*,*) "HIER!! hdmcon()", k
        stop
        hdmcon = 0.0
      end select
C
      return
      end
C
C     *****************************************
C
      subroutine setiprob(N)
C
      integer N
C
      INCLUDE 'commons/odeprb.cmn'
C
      iprob = N
C
      return
      end
C
C     *****************************************
C
      subroutine dfill(N,value,start,step)
C
      integer N, step
      double precision value, start(N)
C
      do ii = 1,N,step
        start(ii) = value
      enddo
C
      return
      end
C
C     *****************************************
C
C     Matrix * Vektor
      subroutine hdmvps(keineAhnung,M,N,MATRIX,LDmatrix,VECTOR,
     &                  OUTPUT,IFERR)
C
      integer N, M, LDmatrix, IFERR, keineAhnung
      double precision MATRIX(M,N), VECTOR(N), OUTPUT(M)
C
      if (keineAhnung .eq. 2) then
        call DGEMV('N',M,N,1.D0,MATRIX,LDmatrix,VECTOR,1,1.D0,OUTPUT,1)
      else 
        call DGEMV('N',M,N,1.D0,MATRIX,LDmatrix,VECTOR,1,0.D0,OUTPUT,1)
      endif
C
c      write(unit=0,fmt=*) MATRIX
c      write(unit=0,fmt=*) VECTOR
c      write(unit=0,fmt=*) OUTPUT
c      write(unit=0,fmt=*) 
C
      return
      end

