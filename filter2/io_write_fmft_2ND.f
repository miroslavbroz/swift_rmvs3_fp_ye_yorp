c  io_write_fmft_2ND.f
c  Write resulting proper FMFT frequencies (for all TPs) to `fmft.out' file.
c  Miroslav Broz (miroslav.broz@email.cz), Apr 8th 2003

      subroutine io_write_fmft_2ND(t,nbod,ntp,istat,g,s,iu,fopenstat)

      include '../swift.inc'
      include 'filter.inc'
      include 'proper_2ND.inc'

      integer iu,nbod,ntp
      integer istat(NTPMAX,NSTAT)
      real*8 t
      real*8 g(-PMAXNPL:PMAXNTP), s(-PMAXNPL:PMAXNTP)
      character*(*) fopenstat

      integer i,ierr,id

      integer i1st
      save i1st
      data i1st /0/

c=======================================================================
c  ... executable code

      if (i1st.eq.0) then
        call io_open(iu,'fmft.out',fopenstat,'FORMATTED',ierr)
        if (ierr.ne.0) then
          write(*,*) ' SWIFT ERROR: in io_write_fmft_2ND: '
          write(*,*) '     Could not open fmft.out:'
          call util_exit(1)
        endif
        i1st = 1
      else
        call io_open(iu,'fmft.out','append','FORMATTED',ierr)
      endif

      do i = 2, nbod
        id = -i
        write(iu,10) id, t/365.25d6, g(id), s(id)
10      format(1x,i5,2x,f23.16,2(2x,1f22.16))
      enddo
 
      do i = 1, ntp
        if (istat(i,1).eq.0) then
          write(iu,10) i, t/365.25d6, g(i), s(i)
        endif
      enddo
 
      close(iu)
 
      return
      end

