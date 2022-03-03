c**********************************************************************
      program makeinpl
c**********************************************************************
c
c  makein creates input file pl.in for numerical integrators,
c  calculates positions and velocities of planets for given elements
c
c  files needed for compilation of this code: 
c  a.f, pleph.f, selcon.f, JPLEPH
c
c  miroslav.broz@email.cz, Mar 13th 2003
c
      implicit none
c  constants
      double precision rad,deg,pi,pi2,eps
c  angles
      parameter(pi = 3.14159265358979d0,
     &  pi2 = 2.d0*pi,
     &  rad = 180.d0/pi,
     &  deg = pi/180.d0)

c  variables read from standard input
      double precision t0,mass,elmts(6)
      integer ncent,m
      logical baryc
      character*255 fmt

c  temporal variables
      integer i,j
      double precision t00,y(6),yp(6),m0,n0

c**********************************************************************
c
c  read input variables
c
c  time
      read(*,*) t0
10    format(a)
c  format of output double precision numbers
      read(*,10) fmt
c  number of planetary bodies
      read(*,*) m

c  equinox 2000.0
      t00=2451545. 
c
c  read masses and elements of m bodies and compute pv
c
      write(*,*) m
      do i=1,m
        read(*,*) mass
        read(*,*) elmts
        do j=3,6
          elmts(j)=elmts(j)*deg
        enddo

c isn't it the Sun in the heliocentric frame?
        if (elmts(1).gt.1.e-16) then
          call pv1(elmts,y(1),y(4),mass)
        else
          do j=1,6
            y(j)=0.d0
          enddo
        endif
        write(*,fmt) mass
        write(*,fmt) y(1),y(2),y(3)
        write(*,fmt) y(4),y(5),y(6)
      enddo

      end

c**********************************************************************
      subroutine pv1(elmts,x,v,gm1)
c**********************************************************************
c
c  calculation of position and velocity from elements
c  elmts(6)   elements of orbit
c  x(3),v(3)  position and velocity in rectangular coords
c  gm1        mass of the planet (or zero)
c
      implicit none
      double precision elmts(6),x(3),v(3),gm1
      double precision a,e,i0,omega,omega0,m0,e0,
     &px,py,pz,qx,qy,qz,ax,ay,az,bx,by,bz,r,gms,kg,gms2,kg2,ekepl1
c      parameter(gms=2.9591220663860403d-4,kg=0.017202098902128322400d0)
      parameter(gms=0.2959122082855911d-03)     ! a new value from JPLEPH

      gms2=gms+gm1
      kg2=sqrt(gms2)

      a=elmts(1)
      e=elmts(2)
      i0=elmts(3)
      omega0=elmts(4)
      omega=elmts(5)
      m0=elmts(6)

      px=cos(omega0)*cos(omega)-sin(omega0)*sin(omega)*cos(i0)
      py=cos(omega0)*sin(omega)+sin(omega0)*cos(omega)*cos(i0)
      pz=sin(omega0)*sin(i0)
      qx=-sin(omega0)*cos(omega)-cos(omega0)*sin(omega)*cos(i0)
      qy=-sin(omega0)*sin(omega)+cos(omega0)*cos(omega)*cos(i0)
      qz=cos(omega0)*sin(i0)

      ax=a*px
      ay=a*py
      az=a*pz
      bx=a*sqrt(1-e**2)*qx
      by=a*sqrt(1-e**2)*qy
      bz=a*sqrt(1-e**2)*qz

      e0=ekepl1(m0,e)

      x(1)=ax*(cos(e0)-e)+bx*sin(e0)
      x(2)=ay*(cos(e0)-e)+by*sin(e0)
      x(3)=az*(cos(e0)-e)+bz*sin(e0)
      r=sqrt(x(1)**2+x(2)**2+x(3)**2)
      v(1)=kg2/(r*sqrt(a))*(-ax*sin(e0)+bx*cos(e0))
      v(2)=kg2/(r*sqrt(a))*(-ay*sin(e0)+by*cos(e0))
      v(3)=kg2/(r*sqrt(a))*(-az*sin(e0)+bz*cos(e0))

      return
      end


