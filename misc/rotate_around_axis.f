c rotate_around_axis.f
c Rotate a vector x around the axis v by the angle phi.
c Miroslav Broz (miroslav.broz@email.cz), Jul 31st 2007

      subroutine rotate_around_axis(v,phi,x,y)

      implicit none

      real*8 v(3), phi, x(3), y(3)
      real*8 c,s

      s = sin(phi)
      c = cos(phi)

      y(1) = x(1) * (c + (1.d0-c)*v(1)**2)       
     :     + x(2) * ((1.d0-c)*v(1)*v(2) - s*v(3))
     :     + x(3) * ((1.d0-c)*v(1)*v(3) + s*v(2))
      y(2) = x(1) * ((1.d0-c)*v(2)*v(1) + s*v(3))
     :     + x(2) * (c + (1.d0-c)*v(2)**2)
     :     + x(3) * ((1.d0-c)*v(2)*v(3) - s*v(1))
      y(3) = x(1) * ((1.d0-c)*v(3)*v(1) - s*v(2))
     :     + x(2) * ((1.d0-c)*v(3)*v(2) + s*v(1))
     :     + x(3) * (c + (1.d0-c)*v(3)**2)

      return
      end


