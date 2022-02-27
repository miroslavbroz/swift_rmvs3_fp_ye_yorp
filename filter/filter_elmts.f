c**********************************************************************
      subroutine filter_elmts(elmts,nbod,ntp,istat,FN,n,filter_M,d)
c**********************************************************************
c
c  filter the array elmts with filter d and store the filtered
c  orbital elements in the level n+1 in array elmts
c
      include '../swift.inc'
      include 'filter.inc'

      real*8 elmts(FMAXN,FMAXE,-FMAXNPL:FMAXNTP,FMAXC)
c  Mar 15th 2002: added istat parameter
      integer nbod,ntp,istat(NTPMAX,NSTAT),n,FN(FMAXC),filter_M
      real*8 d(-FMAXM:FMAXM)

c  functions
      real*8 filter_filter

      integer i,k,id

!$OMP PARALLEL
!$OMP& SHARED (ntp, elmts, FN, n, d, filter_M)
!$OMP DO
!$OMP& PRIVATE (i, k, id)
      do i=2,nbod
        id=-i
        do k=1,FMAXE    ! I DO need sigma for PLs in case of SyMBA
          elmts(FN(n+1),k,id,n+1)=
     &      filter_filter(d,elmts(1,k,id,n),filter_M)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO
!$OMP& PRIVATE (i, k, id)
      do i=1,ntp
        if (istat(i,1).eq.0) then
          id=i
          do k=1,FMAXE
            elmts(FN(n+1),k,id,n+1)=
     &        filter_filter(d,elmts(1,k,id,n),filter_M)
          enddo
        endif
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end       ! filter_elmts


