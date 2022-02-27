
c****************************************************************
      subroutine read_s(iu, s, i)
c****************************************************************
      implicit none
      integer iu,i
      character*(*) s
c  functions
      integer length

5     continue
        i=i+1
        read(iu, 10, end=20, err=20) s
10      format(a)
      if ((s(1:1).eq.'#').or.(length(s).eq.0)) goto 5
      return
20    continue
      write(*,*) 'read_s: Error reading input at line ', i, '.'
      stop
      end

