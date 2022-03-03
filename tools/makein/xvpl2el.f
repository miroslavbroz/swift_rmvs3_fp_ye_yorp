c xvpl2el.f
c Convert positions and velocities in the PL.IN file to orbital elements.
c Miroslav Broz (miroslav.broz@email.cz), Jul 21st 2006

      program xvpl2el

      implicit none
      integer ntp
      real*8 y(6), elmts(6), mass
      real rstat
      integer istat
      integer i, j
      real*8 nula2pi
      real*8 a,h,k,p,q,lambda
      real*8 pi,deg,rad
      parameter(pi = 3.1415926535897932d0, deg = pi/180., rad = 180./pi)

      read(*,*) ntp
      write(*,*) ntp

      do i = 1, ntp

        read(*,*) mass
        read(*,*) y(1), y(2), y(3)
        read(*,*) y(4), y(5), y(6)
c        write(*,*) y(1), y(2), y(3)
c        write(*,*) y(4), y(5), y(6)

        call e1_new(y(1),y(4),elmts,mass)

        write(*,10) elmts(1), elmts(2), (nula2pi(elmts(j))*rad, j=3,6)
10      format(6(1x,f22.16))

c        a = elmts(1)
c        h = elmts(2) * sin(elmts(4)+elmts(5))
c        k = elmts(2) * cos(elmts(4)+elmts(5))
c        p = sin(elmts(3)/2) * sin(elmts(5))
c        q = sin(elmts(3)/2) * cos(elmts(5))
c        lambda = elmts(6) + elmts(4) + elmts(5)

c        write(*,20) a,h,k,p,q,nula2pi(lambda)
c20      format(1x,f12.9,4(1x,f12.9),1x,f12.8)

      enddo

      end

c**********************************************************************
      subroutine e1_new(x,v,elmts,gm1)
c**********************************************************************
c
c  determination of elemnts from the position and velocity
c  x(3)      rectagular coordinates of the body [au]
c  v(3)      velocity [au/day] 
c  elmts(6)  elements of orbit
c  gm1       mass of the planet (or zero)
c
      implicit none
      double precision x(3),v(3),elmts(6),gm1
      double precision v2,r,a,xx,yy,zz,p,e,i0,e0,m0,v0,
     &omega,omega0,sinv,cosv,nula2pi,gms,kg,gms2
c      parameter(gms=2.9591220663860403d-4,kg=0.017202098902128322400d0)
      parameter(gms=0.2959122082855911d-03)     ! a new value from JPLEPH

      gms2=gms+gm1
      v2=v(1)**2+v(2)**2+v(3)**2
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
      a=1/(2/r-v2/gms2)
      xx=x(2)*v(3)-x(3)*v(2)
      yy=x(1)*v(3)-x(3)*v(1)
      zz=x(1)*v(2)-x(2)*v(1)
      p=(xx**2+yy**2+zz**2)/gms2
      e=sqrt(1-p/a)
      i0=atan(sqrt(xx**2+yy**2)/zz)
      omega=nula2pi(atan2(xx,yy))
      sinv=(x(1)*v(1)+x(2)*v(2)+x(3)*v(3))/r*sqrt(p/gms2)
      cosv=(p/r-1)
      v0=nula2pi(atan2(sinv,cosv))
      e0=2*atan(sqrt((1-e)/(1+e))*tan(v0/2))
      m0=nula2pi(e0-e*sin(e0))
      omega0=nula2pi(atan2(x(3)/sin(i0),x(1)*cos(omega)+
     &x(2)*sin(omega))-v0)

      elmts(1)=a
      elmts(2)=e
      elmts(3)=i0
      elmts(4)=omega0
      elmts(5)=omega
      elmts(6)=m0
      return
      end

