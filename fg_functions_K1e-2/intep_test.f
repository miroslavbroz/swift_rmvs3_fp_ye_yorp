      program intep_test

      implicit none

      integer MAX
      parameter(MAX=1000)

      real*8 x(MAX), y(MAX)
      real*8 x0, y0, x1, x2, dx
      integer i, ierr, n, m
      character*255 filename, str

      call getarg(1,filename)
      call getarg(2,str)
      read(str,*) x1
      call getarg(3,str)
      read(str,*) x2
      call getarg(4,str)
      read(str,*) dx

      open(unit=10, file=filename, status='old', iostat=ierr)
      if (ierr.ne.0) then
        write(*,*) "Error opening file '", trim(filename), "'."
        stop
      endif

      i = 0
      ierr = 0
      do while (ierr.eq.0)
        read(10,'(a)',iostat=ierr) str
        if (ierr.ne.0) exit
        if ((str(2:3).ne.'#').and.(len(trim(str)).gt.0).and.(i.lt.MAX))
     :    then
          i = i+1
          read(str,*,iostat=ierr) x(i), y(i)
c          write(*,*) i, x(i), y(i)  ! dbg
        endif
      enddo
      n = i

      x0 = x1
      do while (x0.lt.x2+1.d-8)
        call intep(x0, y0, x, y, n, ierr)
        write(*,*) x0, y0
        x0 = x0+dx
      enddo

      stop
      end


