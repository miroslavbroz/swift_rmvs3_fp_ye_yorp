c vector_product.f
c Vector product of two vectors.
c Miroslav Broz (miroslav.broz@email.cz), Jul 31st 2007

      subroutine vector_product(a,b,c)

      implicit none

      real*8 a(3), b(3), c(3)

      c(1) = a(2)*b(3) - a(3)*b(2);
      c(2) = a(3)*b(1) - a(1)*b(3);
      c(3) = a(1)*b(2) - a(2)*b(1);

      return
      end


