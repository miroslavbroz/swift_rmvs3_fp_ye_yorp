c lsm.f
c Least Square Method, linear fit of y = a*x + b .
c Miroslav Broz (miroslav.broz@email.cz), Aug 18th 2008

      subroutine lsm(x, y, n, a, b, R)

      implicit none

      integer n
      real*8 x(n), y(n), a, b, R

      integer i
      real*8 sx, sxx, sy, syy, sxy,aa,bb
      
      if (n.lt.2) then
        a = 0.d0
        b = 0.d0
        R = 0.d0
        return
      endif

      sx  = 0.d0
      sxx = 0.d0
      sy  = 0.d0
      syy = 0.d0
      sxy = 0.d0

      do i = 1, n
        sx  = sx  + x(i)
        sxx = sxx + x(i)*x(i)
        sy  = sy  + y(i)
        syy = syy + y(i)*y(i)
        sxy = sxy + x(i)*y(i)
      enddo

      a = (n*sxy-sx*sy)/(n*sxx-sx*sx)
      b = (sy-aa*sx)/n
      R = (n*sxy-sx*sy)/sqrt((n*sxx-sx*sx)*(n*syy-sy*sy))

c      a = sy/sx        ! this is for b = 0
c      b = 0.d0

      return
      end


