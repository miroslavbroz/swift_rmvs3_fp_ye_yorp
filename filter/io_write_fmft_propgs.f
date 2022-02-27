c  io_write_fmft_propgs.f
c  Write ALL FMFT frequencies, amplitudes and phases (for all PLs and TPs) to `propgs.out' file.
c  Miroslav Broz (miroslav.broz@email.cz), Nov 30th 2004

      subroutine io_write_fmft_propgs(t,nbod,ntp,istat,
     :  goutput,soutput,iu,fopenstat)

      include '../swift.inc'
      include 'filter.inc'
      include 'proper.inc'

      integer iu,nbod,ntp
      integer istat(NTPMAX,NSTAT)
      real*8 t
      real*8 goutput(-PMAXNPL:PMAXNTP,MAX_FMFT_NFREQ,3)
      real*8 soutput(-PMAXNPL:PMAXNTP,MAX_FMFT_NFREQ,3)
      character*(*) fopenstat

      integer i,j,ierr,id

      integer i1st
      save i1st
      data i1st /0/

c=======================================================================
c  ... executable code

      if (i1st.eq.0) then
        call io_open(iu,'propgs.out',fopenstat,'FORMATTED',ierr)
        if (ierr.ne.0) then
          write(*,*) ' SWIFT ERROR: in io_write_fmft: '
          write(*,*) '     Could not open fmft.out:'
          call util_exit(1)
        endif
        i1st = 1
      else
        call io_open(iu,'propgs.out','append','FORMATTED',ierr)
      endif

      do i = 2, nbod
        id = -i
        do j = 1, fmft_nfreq
          write(iu,10) id, t/365.25d6, j,
     :      goutput(id,j,1), goutput(id,j,2), goutput(id,j,3)*DEGRAD,
     :      soutput(id,j,1),
     :      sin(2.d0*asin(soutput(id,j,2))),
     :      soutput(id,j,3)*DEGRAD
10        format(i5,1x,f23.16,1x,i5,6(1x,f22.16))
        enddo
      enddo
 
      do i = 1, ntp
        if (istat(i,1).eq.0) then
          do j = 1, fmft_nfreq
            write(iu,10) i, t/365.25d6, j,
     :        goutput(i,j,1), goutput(i,j,2), goutput(i,j,3)*DEGRAD,
     :        soutput(i,j,1),
     :        sin(2.d0*asin(soutput(i,j,2))),
     :        soutput(i,j,3)*DEGRAD
          enddo
        endif
      enddo
 
      close(iu)
 
      return
      end


