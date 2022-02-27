
c****************************************************************
      real*8 function arr_avg(a, n)
c****************************************************************
      implicit none
      integer n
      real*8 a(n)

      integer i
      real*8 suma

      suma = 0.d0
      do i = 1, n
        suma = suma + a(i)
      enddo
      arr_avg = suma / n

      return
      end

