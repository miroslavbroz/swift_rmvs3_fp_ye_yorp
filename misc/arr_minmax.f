
c****************************************************************
      real*8 function arr_minmax(a, n, minmax)
c****************************************************************
      implicit none
      integer n
      real*8 a(n)
      logical minmax

      integer i
c  functions
      integer arr_minind, arr_maxind

      if (minmax) then
        i = arr_minind(a, n)
      else
        i = arr_maxind(a, n)
      endif
      arr_minmax = a(i)

      return
      end

