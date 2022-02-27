
c****************************************************************
      integer function read_i(iu,i)
c****************************************************************
      implicit none
      integer iu,i
      character*255 s

      call read_s(iu, s, i)
      read(s, *, end=10, err=10) read_i
      return
10    continue
      write(*,*) 'read_i: Error reading input at line ', i, '.'
      stop
      end

