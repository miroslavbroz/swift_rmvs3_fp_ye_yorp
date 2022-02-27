c****************************************************************
       double precision function filter_filter(d,w,Mm)
c****************************************************************
c
c  filter the data w with filter d
c  d    filter
c  Mn   dimension of the filter and data
c  w    input data
c
      include '../swift.inc'
      include 'filter.inc'
      double precision d(-FMAXM:FMAXM),w(0:FMAXN-1)
      integer Mm

      integer m
      double precision w2

      w2=0.0d0
      do m=-Mm,Mm
        w2=w2+d(m)*w(Mm-m)
      enddo
      filter_filter=w2
      return
      end

