c  baryc2.f
c
c  Apply barycenric correction for indirect effect of given planets
c  (w. ID NUMBERS given on command-line) to the initial conditions
c  of outer ones and asteroids as printed by makein program.
c
c  miroslav.broz@email.cz, Sep 5th 2004
c
      program baryc2
      implicit none

      integer nplmax
      parameter(nplmax = 100)

      integer npl,ntp,nids,id
      integer i,j
      real*8 mass(1:nplmax)
      real*8 pvh(6,1:nplmax)
      real*8 pvht(6)
      integer plids(1:nplmax)
      integer istat
      real*8 rstat
      real*8 massb
      real*8 pvb(6)
      character*80 fmt,str

      integer iargc

      data massb, pvb /0.0d0, 6*0.0d0/
      data fmt /'(3(1x,e22.16))'/
c-----------------------------------------------------------------------
c
c  read planet id's to subtract from the command-line
c
      plids(1)=1
      do i=2,nplmax
        plids(i)=-1
      enddo
      if (iargc().gt.0) then
        nids=min(iargc(),nplmax)
        do i=1,nids
          call getarg(i,str)
          read(str,*) id
          plids(id)=1
        enddo
      else
c subtract inner planets by default
        nids=4
        plids(2)=1
        plids(3)=1
        plids(4)=1
        plids(5)=1
      endif

c-----------------------------------------------------------------------
c
c  read makein output (pl.in and tp.in) from stadard input
c

c  skip comments
5     continue
        read(*,10) str
10      format(a)
      if (str(1:1).eq.'#') goto 5 
      read(str,*) npl

      do i = 1, npl
        read(*,*) mass(i)
        read(*,*) (pvh(j,i), j = 1, 3)
        read(*,*) pvh(4,i),pvh(5,i),pvh(6,i)
      enddo

c
c  calculate the barycenter of the inner solar system
c
      do i = 1, npl
        if (plids(i).gt.0) then
          do j = 1, 6
            pvb(j) = pvb(j) + mass(i) * pvh(j,i)
          enddo
          massb = massb + mass(i)
        endif
      enddo
      do j = 1, 6
        pvb(j) = pvb(j) / massb
      enddo

c      write(*,fmt) (pvb(j), j = 1, 3)
c      write(*,fmt) (pvb(j)*365.25, j = 4, 6)
c      do i = 2, 5
c        write(*,'(f20.10)') mass(1)/mass(i)
c      enddo
c      stop

c
c  apply the barycentric correction to the rest of planets and
c  write everything to the standard output in the same format as input
c
      write(*,*) npl - nids
      write(*,fmt) massb
      write(*,fmt) (pvh(j,1), j = 1, 3)
      write(*,fmt) (pvh(j,1), j = 4, 6)
      do i = 1, npl
        if (plids(i).lt.0) then
          write(*,fmt) mass(i)
          write(*,fmt) (pvh(j,i) - pvb(j), j = 1, 3)
          write(*,fmt) (pvh(j,i) - pvb(j), j = 4, 6)
        endif
      enddo
c
c  do the same for the test particles (do NOT use a big array,
c  because sometimes we need to process ~100000 TPs)
c
      read(*,*) ntp
      write(*,*) ntp                   

      do i = 1, ntp
        read(*,*) pvht(1),pvht(2),pvht(3)
        read(*,*) pvht(4),pvht(5),pvht(6)
        read(*,*) istat
        read(*,*) rstat

        write(*,fmt) (pvht(j) - pvb(j), j = 1, 3)
        write(*,fmt) (pvht(j) - pvb(j), j = 4, 6)
        write(*,*) '0'
        write(*,*) '0.0'
      enddo

      end

