c*************************************************************************
c                            IO_READ_LINE_R8
c*************************************************************************
c read one line from real*4 binary file.
c
c      Input:
c            iu       ==> unit number to write to
c      Output:
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c       Returns:
c      io_read_line_r    ==>   =0 read ok
c                           !=0 read failed is set to iostat variable
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 

      integer function io_read_line_r8(iu,id,a,e,inc,capom,omega,capm) 

      include 'swift.inc'
c      include 'io.inc'

c...  Inputs: 
      integer iu

c...  Output: 
      integer id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      integer*2 id2
      real*8 a8,e8,inc8,capom8,omega8,capm8
      integer ierr

c----
c...  Executable code 

      read(iu,iostat=ierr) id2,a8,e8,inc8,capom8,omega8,capm8
      io_read_line_r8 = ierr
      if(ierr.ne.0) then
         return
      endif

      id = id2

      a = a8
      e = e8
      inc = inc8
      capom = capom8
      capm = capm8
      omega = omega8

      return
      end      ! io_read_line_r
c--------------------------------------------------------------------------

