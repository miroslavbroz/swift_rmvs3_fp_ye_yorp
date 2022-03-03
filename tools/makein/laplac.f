c  Compute total angular momentum of planets from makein.out file
c  and transform coordinates to an invariant plane (ie. Lx = Ly = 0).
c  Miroslav Broz (miroslav.broz@email.cz), Mar 13th 2003

      program laplac

      include 'laplac.inc'

      integer npl, ntp
      real*8 mass(NPLMAX), rh(3,NPLMAX), vh(3,NPLMAX)
      real*8 rht(3), vht(3)

      real*8 Lh(3), Ls(3), x1(3), x2(3)
      character*255 istatstr, rstatstr

      integer i, k, iua, iub
      character*80 fmt, str

c  functions
      integer length

      data fmt /'(6(1x,e22.16))'/

c=======================================================================
c
c  read standard input
c

c  skip comments
5     continue
        read(*,10) str
10      format(a)
      if (str(1:1).eq.'#') goto 5 
      read(str,*) npl

      if (npl.gt.NPLMAX) then
        write(*,20) npl, NPLMAX
20      format('laplac: Number of PLs ', i5, ' .gt. NPLMAX = ', i5, '.')
        stop
      endif
      do i = 1, npl
        read(*,*,err=990,end=990) mass(i)
        read(*,*,err=990,end=990) rh(1,i), rh(2,i), rh(3,i)
        read(*,*,err=990,end=990) vh(1,i), vh(2,i), vh(3,i)
      enddo
c
c  calculate angular momentum and rotation angles
c
      call ang_mom(npl, mass, rh, vh, Lh)

      write(*,'(a,/,a,'//fmt(2:))
     :  '# laplac: Total angular momentum in ecliptic frame '
     :  // 'L_ecl (x, y, z):',
     :  '# ',Lh(1), Lh(2), Lh(3) 

      call cart_sphe(Lh, Ls)
      write(*,'(a,/,a,'//fmt(2:))
     :  '# laplac: T. a. m. in spherical ecliptic coordinates '
     :  // 'L_ecl (r, theta/deg, phi/deg):',
     :  '# ',Ls(1), Ls(2)*rad, Ls(3)*rad 

      call rotat(Lh, Ls(2), Ls(3))
c      write(*,*) '... rotated around z, y axes ', Lh(1), Lh(2), Lh(3)

c
c  transform r and v vectors
c
      do i = 1, npl
        call rotat(rh(1,i), Ls(2), Ls(3))
        call rotat(vh(1,i), Ls(2), Ls(3))
      enddo

      call ang_mom(npl, mass, rh, vh, Lh)

      write(*,'(a,/,a,'//fmt(2:))
     :  '# laplac: Total angular momentum in invariant frame '
     :  // 'L_inv (x, y, z):',
     :  '# ',Lh(1), Lh(2), Lh(3) 
c
c  write standard output for PLs
c
      write(*,*) npl
      do i = 1, npl
        write(*,fmt) mass(i)
        write(*,fmt) rh(1,i), rh(2,i), rh(3,i)
        write(*,fmt) vh(1,i), vh(2,i), vh(3,i)
      enddo
c
c  read and write TPs
c
      read(*,*) ntp
      write(*,*) ntp

      do i = 1, ntp
        read(*,*) rht(1), rht(2), rht(3)
        read(*,*) vht(1), vht(2), vht(3)
        read(*,*) istatstr
        read(*,*) rstatstr

c  do the transformation to the invariant plane
        call rotat(rht, Ls(2), Ls(3))
        call rotat(vht, Ls(2), Ls(3))

        write(*,fmt) rht(1), rht(2), rht(3)
        write(*,fmt) vht(1), vht(2), vht(3)
        write(*,10) istatstr(1:length(istatstr))
        write(*,10) rstatstr(1:length(rstatstr))
      enddo

      stop

990   continue
      write(*,*) 'laplac: Stdin format is broken.'

      end

c=======================================================================
c  Vectorial product c = a x b.

      subroutine vdot(a, b, c)
      real*8 a(3), b(3), c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
      return
      end

c=======================================================================
c  Total angular momentum L = Sum_i^n r_i x m_i v_i.

      subroutine ang_mom(n, m, r, v, L)

      include 'laplac.inc'

      integer n
      real*8 m(NPLMAX), r(3,NPLMAX), v(3,NPLMAX), L(3)

      integer i, k
      real*8 x1(3), x2(3)

      do k = 1, 3
        L(k) = 0.d0
      enddo
      do i = 1, n
        do k = 1, 3
          x1(k) = m(i) * v(k, i)
        enddo
        call vdot(r(1,i), x1, x2)
        do k = 1, 3
          L(k) = L(k) + x2(k)
        enddo
      enddo

      return
      end

c=======================================================================
c  Transformation: cartesian -> spherical coordinates.
c  b(1) = r, b(2) = phi, b(3) = theta in usual notation.

      subroutine cart_sphe(a, b)

      implicit none
      real*8 a(3), b(3)

      b(1) = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      b(2) = datan2(a(2), a(1))
      b(3) = dasin(a(3)/b(1))

      return
      end

c=======================================================================
c  Rotation around z and y axes.

      subroutine rotat(x, a, b)

      include 'laplac.inc'

      real*8 x(3), a, b

      real*8 x0(3)

      call rotat3(x, x0, a)
      call rotat2(x0, x, 0.5d0*pi - b)

      return
      end

