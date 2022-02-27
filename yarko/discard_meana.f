c**********************************************************************
c DISCARD_MEANA.F
c**********************************************************************
c Discard TP's according to their MEAN semimajor axis (depends
c on filtering, see filter.in file)
c
c Input:
c  t			time
c  ntp			number of test particles
c  mean_a1, mean_a2	criteria for mean semimajor axis [AU]
c  istat(NTPMAX,NSTAT)	flag in istat array (istat(i,2).eq.-6)
c  rstat(NTPMAX,NSTAT)	if discarded, rstat(i,2) = mean a
c  mean orbital elements in common block /mean_elements/
c
c Remarks: 
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Nov 15th 2000

      subroutine discard_meana(t,ntp,a_mean1,a_mean2,istat,rstat)

      include '../swift.inc'
      include '../filter/filter.inc'

c  input/output
      integer ntp,istat(NTPMAX,NSTAT)
      real*8 t,a_mean1,a_mean2,rstat(NTPMAX,NSTAT)

c  input orbital elements last computed in filter_elmts_write.f subroutine
c  in common block /mean_elements/, we do NOT actually need the mean anomaly
c  (nor peri or node); the upflg checks, if elements were just updated
      real*8 mean_elmts(3,1:FMAXNTP)
      logical upflg
      common/mean_elements/mean_elmts,upflg

c  internal
      integer i
      real*8 a
      integer i1st
      data i1st /1/

c  main ----------------------------------------------------------
      if ((a_mean1.gt.0).or.(a_mean2.gt.0)) then
        if (i1st.eq.1) then
          upflg = .false.
          i1st = 0
        endif
        if (.not.upflg) return	! exit unless mean elements were updated

        do i=1,ntp
          if (istat(i,1).eq.0) then
            a = mean_elmts(1,i)
            if (((a.lt.a_mean1).and.(a_mean1.gt.0)).or.
     :        ((a.gt.a_mean2).and.(a_mean2.gt.0))) then
              write(*,*) 'Particle',i,' mean a=',a,' AU,',
     :          ' discarded at t=',t
              istat(i,1) = 1
              istat(i,2) = -6
              rstat(i,2) = a
            endif
          endif
        enddo

        upflg = .false.		! TP's were discarded, wait for next update
      endif

      return
      end	! discard_meana.f
c-----------------------------------------------------------------

