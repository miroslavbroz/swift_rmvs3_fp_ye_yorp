
c****************************************************************
      logical function ang_eq(a, b, toler)
c****************************************************************
c
c  Test, if two angles <0, 2pi) are equal within tolerance.
c
      real*8 a, b, toler

      real*8 PI,TWOPI
      parameter (PI = 3.14159265358979D0)
      parameter (TWOPI = 2.0D0 * PI)

      if ((abs(a-b).lt.toler).or.
     :  (abs(a-b-TWOPI).lt.toler).or.
     :  (abs(a-b+TWOPI).lt.toler)) then
        ang_eq = .true.
      else
        ang_eq = .false.
      endif

      return
      end

