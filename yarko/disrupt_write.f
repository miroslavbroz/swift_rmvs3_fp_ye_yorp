c**********************************************************************
c DISRUPT_WRITE.F
c**********************************************************************
c Write disrupt.out file with current disruption flags.
c
c Input:
c  t			time
c  ntp			number of test particles
c  istat(NTPMAX,NSTAT)	disruption flag in istat array (istat(i,2).eq.-5)
c  rstat(NTPMAX,NSTAT)	not yet used
c  iu			unit for writing disrupt.out
c  fopenstat		status of disrupt.out
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Nov 8th 2000

      subroutine disrupt_write(outfile,t,ntp,istat,rstat,iu,fopenstat)

      include '../swift.inc'
      integer ntp,istat(NTPMAX,NSTAT),iu
      real*8 t,rstat(NTPMAX,NSTAT)
      character*(*) outfile
      character*80 fopenstat

      integer i,ierr,i1st
      logical disrupt_flag(NTPMAX)
      data i1st/0/
      save i1st

      if (i1st.eq.0) then
        call io_open(iu,outfile,fopenstat,'formatted',ierr)
        if(ierr.ne.0) then
          write(*,*) ' SWIFT ERROR: in disrupt_write'
          write(*,*) '     Could not open ', outfile
          call util_exit(1)
        endif
        i1st = 1
      else
        call io_open(iu,outfile,'append','formatted',ierr)
      endif

      do i=1,ntp
        if ((istat(i,1).eq.1).and.(istat(i,2).eq.-5)) then
          disrupt_flag(i) = .true.
        else
          disrupt_flag(i) = .false.
        endif
      enddo

      write(iu,10) t/365.25d6,(disrupt_flag(i),i=1,ntp)
10    format(1x,1e17.10,1000(1x,l1))

      close(unit = iu)

      return
      end	! disrupt_write.f
c-----------------------------------------------------------------


