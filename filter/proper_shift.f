
c****************************************************************
      subroutine proper_shift(tstart,tstop,PN,nbod,ntp,
     :  prop_time,prop_elmts,nelmts)
c****************************************************************
c
c  shift the prop_* buffers back in time
c
      include '../swift.inc'
      include 'filter.inc'
      include 'proper.inc'

      integer PN,nbod,ntp,nelmts
      real*8 tstart,tstop
      real*8 prop_time(PMAXN)
      real*8 prop_elmts(PMAXN,PMAXE,-PMAXNPL:PMAXNTP)

      integer i,j,k,l,id

      tstop = tstop + prop_dt
      tstart = tstop - prop_win
      if (tstart.ge.prop_time(PN)) then
        PN = 0
      else
        k = 1
        do while (prop_time(k).lt.tstart)
          k = k + 1
        enddo
        call arr_shift(prop_time, PN, k)
        do i = 2, nbod
          id = -i
          do j = 1, nelmts
            call arr_shift(prop_elmts(1,j,id), PN, k)
          enddo 
        enddo
        do i = 1, ntp
          id = i
          do j = 1, nelmts
            call arr_shift(prop_elmts(1,j,id), PN, k)
          enddo 
        enddo
        PN = PN - k
      endif

      return
      end


