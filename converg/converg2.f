c converg.f
c Check convergence of orbits wrt. the median orbit.
c Miroslav Broz (miroslav.broz@email.cz), Jun 26th 2021

      subroutine converg2(time,ntp,mass1,xht,yht,zht,vxht,vyht,vzht,
     :  istat,iu)

      include "../swift.inc"

      real*8 day, au
      parameter(day = 86400.d0)
      parameter(au = 149597870700.d0)

c input
      real*8 time
      real*8 mass1
      integer ntp
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      integer istat(ntp)
      integer iu

c internal
      integer i,j,k,l,n,ialpha
      real*8 a(NTPMAX),e(NTPMAX),inc(NTPMAX),capom(NTPMAX),
     :  omega(NTPMAX),capm(NTPMAX)
      real*8 Dv(NTPMAX),Dvori(NTPMAX),Dvbst(NTPMAX)
      real*8 amean,nmean,DOmega,Dvarpi,val
      real*8 Dvmed,Dvmin,Dvmax
      real*8 a1,capom1,omega1
      integer iub,iuc,ierr
      save iub
      logical regular(NTPMAX)
      save regular

c functions
      real*8 f360,median,median_and_sort
      integer minabsind

      integer i1st
      save i1st
      data i1st /1/

      integer wrt,what
      real*8 vlimit
      real*8 t1,t2
      integer nclones
      save vlimit,t1,t2,nclones

c read parameters
      if (i1st.eq.1) then
        write(*,*) 'Enter with-respect-to body : '
        read(*,*) wrt
        write(*,*) 'Enter what body : '
        read(*,*) what
        write(*,*) 'Enter velocity limit in m/s : '
        read(*,*) vlimit 
        write(*,*) 'Enter time span t1, t2 : '
        read(*,*) t1,t2
        write(*,*) 'Enter nclones : '
        read(*,*) nclones

        write(*,*) 'wrt = ', wrt
        write(*,*) 'what = ', what
        write(*,*) 'vlimit = ', vlimit
        write(*,*) 't1 = ', t1
        write(*,*) 't2 = ', t2
        write(*,*) 'nclones = ', nclones

        open(unit=iu, file="converg2.out", status="unknown")
        write(iu,*) '# t    Dvmed Dvmin Dvmax Dv(', what, ') j '
        write(iu,*) '# [My] [m/s] [m/s] [m/s] [m/s] []'

        iub=105
        open(unit=iub, file="converg2.tmp", status="unknown")
        write(iub,*) '# id t    Dv(i)'
        write(iub,*) '# [] [My] [m/s]'

        iuc=110
        open(unit=iuc, file="regular.lst", status="old")
        do i = 1,ntp
          read(iuc,*,iostat=ierr) j,regular(i)
          if ((ierr.ne.0).or.(j.ne.i)) then
            write(*,*) 'Error reading regular.lst.'
            stop
          endif
          write(*,*) i,regular(i)
        enddo
        close(iuc)

        i1st=0
      endif

      if (istat(1).ne.0) then
        write(*,*) 'converg.f: Warning, test particle 1 not present!'
        return
      endif

      do i = 1,ntp
        if (istat(i).eq.0) then
          call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     &      vzht(i),mass1,ialpha,a(i),e(i),inc(i),capom(i),omega(i),
     &      capm(i))
        endif
      enddo

c 2DO use istat or regular.lst!
c compute median orbit
      a1 = median(ntp,a)
      capom1 = median(ntp,capom)
      omega1 = median(ntp,omega)

c deltav's wrt. median orbit
      do i = 1,ntp
        if ((istat(i).eq.0).and.(regular(i))) then
          amean = 0.5d0*(a(i)+a1)
          nmean = sqrt(mass1/amean**3)
          DOmega = f360(capom(i)-capom1)
          Dvarpi = f360(capom(i)+omega(i)-(capom1+omega1))
          val = nmean*amean*sqrt((sin(inc(i))*DOmega)**2
     :      + (e(i)*Dvarpi)**2)
          Dv(i) = val
          Dvori(i) = val
        endif
      enddo

c choose best clone
      n = ntp/nclones
      do i = 1,n
        j = (i-1)*nclones+1
        l = j-1 + minabsind(nclones,Dv(j))
        Dvbst(i) = Dv(l)
      enddo

c compute statistics (do not use 1, which is always close to median)
      Dvmed = median_and_sort(n, Dvbst)
      Dvmin = Dvbst(2)
      Dvmax = Dvbst(n)

c output w. fine sampling
      if (Dvmed.lt.vlimit/(au/day)) then
        write(iu,*) time/365.25d6, Dvmed*(au/day), Dvmin*(au/day),
     :    Dvmax*(au/day), Dvori(what)*(au/day), n
      endif

      if ((time.gt.t1).and.(time.lt.t2)) then
        do i = 1,j
          write(iub,*) i, time/365.25d6, Dvori(i)*(au/day)
        enddo
        write(iub,*) 
      endif

      return
      end

      integer function minabsind(n, x)
      implicit none
      integer n
      real*8 x(n)
      integer i,j
      j = 1
      do i = 2, n
        if (abs(x(i)).lt.abs(x(j))) then
          j = i
        endif
      enddo
      minabsind = j
      return
      end


