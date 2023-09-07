c dahlgren2.f
c Check close encounters wrt. a planet "wrt".
c Miroslav Broz (miroslav.broz@email.cz), Jun 26th 2021

c Reference: Dahlgren (1998), A&A 336, 1056-1064.

c Note: Collisional probabilities are computed ex-post.

      subroutine dahlgren2(time,ntp,nbod,mass1,xh,yh,zh,vxh,vyh,vzh,
     :  xht,yht,zht,vxht,vyht,vzht,istat,iu)

      include "../swift.inc"

      real*8 day, au
      parameter(day = 86400.d0)
      parameter(au = 149597870700.d0)

c input
      real*8 time
      real*8 mass1
      integer ntp,nbod
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      integer istat(ntp)
      integer iu

c internal
      integer i,j
      real*8 dx,dy,dz,r
      real*8 dvx,dvy,dvz,v

      integer i1st
      save i1st
      data i1st /1/

      real*8 rlimit
      integer wrt
      logical use_dahlgren
      save rlimit,wrt,use_dahlgren

c read parameters
      if (i1st.eq.1) then
        write(*,*) 'rlimit = '
        read(*,*) rlimit
        write(*,*) rlimit

        write(*,*) 'wrt = '
        read(*,*) wrt
        write(*,*) wrt

        write(*,*) 'use_dahlgren  = '
        read(*,*) use_dahlgren
        write(*,*) use_dahlgren

        open(unit=iu, file="dahlgren.out", status="unknown",
     :    access="append")
        write(iu,*) '# t i j r v'
        write(iu,*) '# My - - au m/s'

        i1st=0
      endif

      if (.not.use_dahlgren) return

c check close encouters
      do i = 1,ntp
        if (istat(i).eq.0) then
          j = wrt
          dx = xht(i)-xh(j)
          dy = yht(i)-yh(j)
          dz = zht(i)-zh(j)
          r = sqrt(dx**2 + dy**2 + dz**2)

          if (r.le.rlimit) then
            dvx = vxht(i)-vxh(j)
            dvy = vyht(i)-vyh(j)
            dvz = vzht(i)-vzh(j)
            v = sqrt(dvx**2 + dvy**2 + dvz**2)

            write(iu,*) time/365.25d6, i, j, r, v*(au/day)
          endif
        endif
      enddo

      return
      end


