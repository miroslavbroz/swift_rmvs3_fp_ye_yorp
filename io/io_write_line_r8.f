c*************************************************************************
c                            IO_WRITE_LINE_R8
c*************************************************************************
c write out one line to real*8 binary file.
c
c      Input:
c            iu       ==> unit number to write to
C	     a        ==> semi-major axis or pericentric distance if a parabola
c                          (real scalar)
c            e        ==> eccentricity (real scalar)
C            inc      ==> inclination  (real scalar)
C            capom    ==> longitude of ascending node (real scalar)
C	     omega    ==> argument of perihelion (real scalar)
C	     capm     ==> mean anomoly(real scalar)
c
c Remarks: 
c Authors:  Hal Levison 
c Date:    2/22/94
c Last revision: 
c
c Modified: Miroslav Broz
c Date: May 2nd 1999

      subroutine io_write_line_r8(iu,id,a,e,inc,capom,omega,capm) 

      include '../swift.inc'
c      include 'io.inc'

c...  Inputs: 
      integer iu,id
      real*8 a,e,inc,capom,omega,capm

c...  Internals
      integer*2 id2
      real*8 a8,e8,inc8,capom8,omega8,capm8


c----
c...  Executable code 

      id2 = id

      a8 = a
      e8 = e
      inc8 = inc
      capom8 = capom
      capm8 = capm
      omega8 = omega

      write(iu) id2,a8,e8,inc8,capom8,omega8,capm8

      return
      end      ! io_write_line_r8
c--------------------------------------------------------------------------
