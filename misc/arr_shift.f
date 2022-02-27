
c****************************************************************
      subroutine arr_shift(a, n, k)
c****************************************************************
      implicit none
      integer n, k
      real*8 a(n)

      integer i

      do i = k+1, n
        a(i-k) = a(i)
      enddo

      return
      end

