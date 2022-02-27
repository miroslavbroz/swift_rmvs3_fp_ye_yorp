c normalize.f
c Normalize a 3-dimensional vector.
c Miroslav Broz (miroslav.broz@email.cz), Aug 6th 2007

      subroutine normalize(x)

      implicit none
      real*8 x(3)

      integer i
      real*8 norm

      norm = 0.d0
      do i = 1, 3
        norm = norm + x(i)*x(i)
      enddo
      norm = sqrt(norm)

      if (abs(norm).gt.tiny(norm)) then
        norm = 1.d0/norm
        x = x*norm
      else
        write(*,'(a)') '# Warning: zero vector in normalize'
      endif

      return
      end


