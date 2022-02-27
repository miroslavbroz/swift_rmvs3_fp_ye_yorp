
c****************************************************************
      integer function isgn(x)
c****************************************************************
c
c  Return an integer signum of a double precision argument.
c
      real*8 x
 
      if (x.gt.0.d0) then
        isgn = 1
      elseif (x.lt.0.d0) then
        isgn = -1
      else
        isgn = 0
      endif

      return
      end

c****************************************************************

