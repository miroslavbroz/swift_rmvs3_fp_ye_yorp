c  Apply the inverse rotation from the invarinat plane to ecliptic
c  reference frame to a text file formatted like follow.out, using
c  angle values from laplac.out.
c  Miroslav Broz (miroslav.broz@email.cz), Sep 30th 2003

      program lapinv

      include 'laplac.inc'

      real*8 a,e,inc,capom,omega,capm
      real*8 r(3), v(3), gm, Ls2, Ls3, t
      integer ialpha, id
      character*80 str,fmt

c=======================================================================
c
c  read standard input
c

c output format string
      fmt = '(1x,i4,1x,f16.10,1x,f14.8,1x,f10.8,4(1x,f10.5))'

c  skip comments
5     continue
        read(*,10,end=990,err=990) str
10      format(a)
      if (str(1:1).eq.'#') goto 5 
      read(str,*) gm, Ls2, Ls3
      Ls2=Ls2*deg
      Ls3=Ls3*deg

c read orbital elements, convert them to cartesian coordinates (r, v),
c rotate, convert back and write output
20    continue
        read(*,*,end=990,err=990) id,t,a,e,inc,capom,omega,capm

        inc=inc*deg
        capom=capom*deg
        omega=omega*deg
        capm=capm*deg
        ialpha=-1
        call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     :    r(1),r(2),r(3),v(1),v(2),v(3))

        call rotat_inv(r,Ls2,Ls3)
        call rotat_inv(v,Ls2,Ls3)

        call orbel_xv2el(r(1),r(2),r(3),v(1),v(2),v(3),gm,ialpha,
     :    a,e,inc,capom,omega,capm)

        write(*,fmt) id,t,a,e,inc*rad,capom*rad,omega*rad,capm*rad
      goto 20

990   continue

      stop
      end

c=======================================================================
c  INVERSE rotation around z and y axes.

      subroutine rotat_inv(x, a, b)

      include 'laplac.inc'

      real*8 x(3), a, b

      real*8 x0(3)

      call rotat2(x, x0, -(0.5d0*pi - b))
      call rotat3(x0, x, -a)

      return
      end

