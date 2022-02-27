
c****************************************************************
      logical function read_l(iu,i)
c****************************************************************
      implicit none
      integer iu,i
      character*255 s

      call read_s(iu, s, i)
      read(s, *, end=10, err=10) read_l
      return
10    continue
      write(*,*) 'read_l: Error reading input at line ', i, '.'
      stop
      end

