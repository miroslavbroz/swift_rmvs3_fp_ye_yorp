
c****************************************************************
      integer function arr_maxind(a, n)
c****************************************************************
      implicit none
      integer n
      real*8 a(n)

      integer i, maxi
      real*8 maxa

      maxa = a(1)
      maxi = 1
      do i = 2, n
        if (a(i).gt.maxa) then
          maxa = a(i)
          maxi = i
        endif
      enddo
      arr_maxind = maxi

      return
      end

