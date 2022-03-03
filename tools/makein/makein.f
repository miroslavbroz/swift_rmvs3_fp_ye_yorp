c**********************************************************************
      program makein
c**********************************************************************
c
c  makein creates input files for numerical integrators,
c  prints positions and velocities of planets in given time and
c  calculates the same for zero mass bodies for given elements
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

c  Mar 13th 2003: commented inclination of ecliptic
c      parameter(eps=0.409099609d0)

c  variables read from standard input
      double precision t0,elmts(6)
      integer ncent,m
      logical baryc
      character*255 fmt

c  variables read from JPLEPH file
      common/masses/gm
      double precision gm(13),gms,kg
      integer index(11)
      character*6 nams(12)

c  temporal variables
      integer i,j
      double precision t00,y(6),yp(6),m0,n0

c  functions
      double precision eps_earth

c  variables read from JPLEPH file
      data index/11,1,2,13,4,5,6,7,8,9,10/
      data nams/'GM1','GM2','GMB','GM4','GM5','GM6','GM7',
     &'GM8','GM9','GM10','GMS','EMRAT'/

c**********************************************************************
c
c  read input variables
c
c  time
      read(*,*) t0
10    format(a)
c  number of central body: 11 Sun, 12 solar system barycenter
      read(*,*) ncent
c  position of Earth/Moon barycenter (T) or Earth alone (F)?
      read(*,*) baryc
c  format of output double precision numbers
      read(*,10) fmt
c  number of zero mass bodies
      read(*,*) m

c  equinox 2000.0
      t00=2451545. 
c
c  estimate inclination of ecliptic to equator
c
      eps = eps_earth((t0-t00)/36525.d0)
c
c  read GM of Sun and planets from JPLEPH file
c
      call selcon(nams,12,gm)
      gms=gm(index(1))
      kg=sqrt(gms)
      write(*,15) kg
15    format('# k_gauss = ', e22.16)

c  Earth/Moon barycenter option
      gm(13)=gm(3)
      if (.not.baryc) then
        index(4)=3
        gm(3)=gm(3)/(1.d0+1.d0/gm(12))
      endif
c
c  read initial positions and velocities of planets from JPLEPH file
c
      write(*,*) '10'
      do i=1,10
        call pleph(t0,index(i),ncent,yp(1))
        call rotat1(yp(1),y(1),-eps)
        call rotat1(yp(4),y(4),-eps)
        write(*,fmt) gm(index(i))
        write(*,fmt) y(1),y(2),y(3)
        write(*,fmt) y(4),y(5),y(6)
      enddo
c
c  read elements of m bodies and compute pv
c
      write(*,*) m
      do i=1,m
        read(*,*) elmts
        do j=3,6
          elmts(j)=elmts(j)*deg
        enddo

c  Mar 10th 1999: commented out J2000.0 equinox
c        m0=elmts(6)
c        n0=kg*elmts(1)**(-1.5)
c        elmts(6)=m0+n0*(t0-t00)

        call pv(elmts,y(1),y(4))
        write(*,fmt) y(1),y(2),y(3)
        write(*,fmt) y(4),y(5),y(6)
        write(*,*) '0'
        write(*,*) '0.0'
      enddo

      end

