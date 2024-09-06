c converg.f
c Check convergence of orbits wrt. particle "wrt".
c Miroslav Broz (miroslav.broz@email.cz), Jun 26th 2021

      subroutine converg(time,ntp,mass1,xht,yht,zht,vxht,vyht,vzht,
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
      integer i,j,k,ialpha
      real*8 a(NTPMAX),e(NTPMAX),inc(NTPMAX),capom(NTPMAX),
     :  omega(NTPMAX),capm(NTPMAX)
      real*8 Dv(NTPMAX),Dvori(NTPMAX)
      real*8 amean,n,DOmega,Dvarpi,val
      real*8 Dvmed,Dvmin,Dvmax
      integer iub,iuc,ierr
      save iub
      logical regular(NTPMAX)
      save regular

c functions
      real*8 f360,median,min2,max2

      integer i1st
      save i1st
      data i1st /1/

      integer wrt,what
      real*8 vlimit,t1,t2
      save wrt,what,vlimit,t1,t2

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

        write(*,*) 'wrt = ', wrt
        write(*,*) 'what = ', what
        write(*,*) 'vlimit = ', vlimit
        write(*,*) 't1 = ', t1
        write(*,*) 't2 = ', t2

        open(unit=iu, file="converg.out", status="unknown")
        write(iu,*) '# t    Dvmed Dvmin Dvmax Dv(', what, ') j '
        write(iu,*) '# [My] [m/s] [m/s] [m/s] [m/s] []'

        iub=105
        open(unit=iub, file="converg.tmp", status="unknown")
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

c check convergence
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

      j = 0
      k = wrt
      do i = 1,ntp
        if ((i.ne.k).and.(istat(i).eq.0).and.(regular(i))) then
          amean = 0.5d0*(a(i)+a(k))
          n = sqrt(mass1/amean**3)
          DOmega = f360(capom(i)-capom(k))
          Dvarpi = f360(capom(i)+omega(i)-(capom(k)+omega(k)))
          val = n*amean*sqrt((sin(inc(i))*DOmega)**2+(e(i)*Dvarpi)**2)
          j = j+1
          Dv(j) = val
          Dvori(j) = val
        endif
      enddo

      Dvmed = median(j, Dv)
      Dvmin = Dv(1)
      Dvmax = Dv(j)

!      if (Dvmed.lt.vlimit/(au/day)) then
      if (Dvmin.lt.vlimit/(au/day)) then
        write(iu,*) time/365.25d6, Dvmed*(au/day), Dvmin*(au/day),
     :    Dvmax*(au/day), Dvori(what)*(au/day), j
      endif

      if ((time.gt.t1).and.(time.lt.t2)) then
        do i = 1,j
          write(iub,*) i, time/365.25d6, Dvori(i)*(au/day)
        enddo
        write(iub,*) 
      endif

      return
      end


      real*8 function f360(x)
      include '../swift.inc'
      real*8 x
      if (x.lt.-PI) then
        x = x + TWOPI
      else if (x.gt.PI) then
        x = x - TWOPI
      endif
      f360 = x
      return
      end

      real*8 function median(n, x)
      implicit none
      integer n
      real*8 x(n)
      real*8 y(n)
      real*8 median_and_sort
      integer i
      do i = 1,n
        y(i) = x(i)
      enddo
      median = median_and_sort(n, y)
      return
      end

      real*8 function median_and_sort(n, x)
      implicit none
      integer n
      real*8 x(n)
      integer j
      real*8 s
      if (n.eq.0) then
        median_and_sort = -1.d0
        return
      endif
      call quicksort(n, x)
      j = int(n/2.d0)
      if (mod(n, 2).eq.0) then
        s = 0.5d0*(x(j)+x(j-1))
      else
        s = x(j)
      endif
      median_and_sort = s
      return
      end


