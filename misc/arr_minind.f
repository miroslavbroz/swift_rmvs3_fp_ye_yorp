
c****************************************************************
      integer function arr_minind(a, n)
c****************************************************************
      implicit none
      integer n
      real*8 a(n)

      integer i, mini
      real*8 mina

      mina = a(1)
      mini = 1
      do i = 2, n
        if (a(i).lt.mina) then
          mina = a(i)
          mini = i
        endif
      enddo
      arr_minind = mini

      return
      end

