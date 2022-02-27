
c****************************************************************
      integer function ang_sgndot(a, b)
c****************************************************************
c
c  Signum of d angle/d t from two consecutive values a, b,
c  asure angles near 0 and 2pi.
c
      real*8 a, b

      real*8 PI,TWOPI
      parameter (PI = 3.14159265358979D0)
      parameter (TWOPI = 2.0D0 * PI)

      integer isgn

      if (abs(b-a-TWOPI).lt.abs(b-a)) then
        ang_sgndot = isgn(b-a-TWOPI)
      else if (abs(b-a+TWOPI).lt.abs(b-a)) then
        ang_sgndot = isgn(b-a+TWOPI)
      else
        ang_sgndot = isgn(b-a)
      endif

      return
      end

