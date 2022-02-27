c**********************************************************************
c REORIENT_WRITE.F
c**********************************************************************
c Write reorient.out file with current spin axes orientations.
c
c  outfile              output filename
c  t			time
c  ntp			number of test particles
c  iu			unit for writing reorient.out
c  fopenstat		status of reorient.out
c  s(3,NTPMAX)		orientations of spin axes (in common block /yarko/)
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Jan 19th 2010

      subroutine reorient_write(outfile,t,ntp,istat,iu,fopenstat)

      include '../swift.inc'
      include 'spin.inc'

      integer ntp,iu,istat(NTPMAX,NSTAT)
      real*8 t
      character*(*) outfile
      character*80 fopenstat

      integer i,ierr,i1st
      data i1st/0/
      save i1st

      if (i1st.eq.0) then
        call io_open(iu,outfile,fopenstat,'formatted',ierr)
        if(ierr.ne.0) then
          write(*,*) ' SWIFT ERROR: in reorient_write'
          write(*,*) '     Could not open ', outfile
          call util_exit(1)
        endif
        i1st = 1
      else
        call io_open(iu,outfile,'append','formatted',ierr)
      endif

      do i = 1,ntp
        if (istat(i,1).eq.0) then
          write(iu,20) i, t/365.25d6, s(1,i),s(2,i),s(3,i), omega(i)
20        format(i5,1x,f18.10,1x,4(1x,1f11.8))
        endif
      enddo

      close(unit = iu)

      return
      end	! reorient_write.f
c-----------------------------------------------------------------


