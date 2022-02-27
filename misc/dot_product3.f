c dot_product3.f
c Dot product of two vectors.
c Miroslav Broz (miroslav.broz@email.cz), Jul 31st 2007

      real*8 function dot_product3(a, b)

      implicit none

      real*8 a(3), b(3)

      dot_product3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

      return
      end


