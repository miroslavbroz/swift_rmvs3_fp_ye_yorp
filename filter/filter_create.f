c****************************************************************
      subroutine filter_create(d,Mm,x0,beta)
c****************************************************************
c
c  create filter
c  d(Mm)       filter
c  Mm,x0,beta  parameters of the filter
c
      include '../swift.inc'
      include 'filter.inc'

c  parameters
      integer Mm
      double precision x0,beta,d(-FMAXM:FMAXM)
c  local varibles
      integer m
      double precision C
c  functions
      double precision bessi0

      d(0)=2.d0*x0*bessi0(beta)
      do m=1,Mm      
        d(m)=sin(2.d0*pi*m*x0)/(pi*m)
     :    * bessi0(beta*sqrt(1.d0-m*m/dble(Mm*Mm)))
        d(-m)=d(m)
      enddo
      C=0.d0
      do m=-Mm,Mm
        C=C+d(m)
      enddo
      C=1.d0/C
      do m=-Mm,Mm
        d(m)=C*d(m)
      enddo
      return
      end

