c**********************************************************************
      subroutine filter_shift(elmts,nbod,ntp,n,FBUFN,filter_D)
c**********************************************************************
c
c  shift left elements in array elmts by filter_D(n) / decimation
c
      include '../swift.inc'
      include 'filter.inc'

      real*8 elmts(FMAXN,FMAXE,-FMAXNPL:FMAXNTP,FMAXC)
      integer nbod,ntp,n,FBUFN,filter_D(FMAXC)

      integer i,j,k,id

!$OMP PARALLEL
!$OMP& SHARED (ntp, filter_D, FBUFN, elmts)
!$OMP DO
!$OMP& PRIVATE (i, j, k, id)
      do i=2,nbod
        id=-i
        do k=1,FMAXE    ! I DO need sigma for PLs in case of SyMBA
          do j=filter_D(n)+1,FBUFN
            elmts(j-filter_D(n),k,id,n)=elmts(j,k,id,n)
          enddo
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO
!$OMP& PRIVATE (i, j, k, id)
      do i = 1, ntp
        id=i
        do k=1,FMAXE
          do j=filter_D(n)+1,FBUFN
            elmts(j-filter_D(n),k,id,n)=elmts(j,k,id,n)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end       ! filter_shift

