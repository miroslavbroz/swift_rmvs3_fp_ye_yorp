c  <0,2pi) truncation

      real*8 function zero2pi(x)
      real*8 x

      real*8 y,PI,TWOPI
      parameter (PI = 3.14159265358979D0, PI2 = 2.0D0 * PI)

      y=x/PI2
      y=(y-aint(y))*PI2
      if (y.lt.0.d0) y=y+PI2
      zero2pi=y

      return
      end

