c**********************************************************************
      subroutine filter_elmts_write(iu,elmts,n,nbod,ntp,istat,d,FM,
     :  filter_direct,filter_write,filter_write_r8)
c**********************************************************************
c
c  filter the array elmts with filter d and write filtered
c  orbital elements to the binary output file iu
c
      include '../swift.inc'
      include 'filter.inc'
      include 'cb_meanel.inc'

      integer iu,nbod,ntp,n,FM,istat(NTPMAX,NSTAT)
      real*8 elmts(FMAXN,FMAXE,-FMAXNPL:FMAXNTP,FMAXC)
      real*8 d(-FMAXM:FMAXM)
      logical filter_direct,filter_write,filter_write_r8
c temp vars
      integer i,id
      real*8 a,k,h,q,p,e,inc,peri,node,varpi,sigma,s_c,s_s
c  functions
      real*8 filter_filter, zero2pi

c  I did'nt parallelize these 2 cycles, because of i/o routines.
      do i=2,nbod
        id=-i
        a = filter_filter(d,elmts(1,1,id,n),FM)
        k = filter_filter(d,elmts(1,2,id,n),FM)
        h = filter_filter(d,elmts(1,3,id,n),FM)
        q = filter_filter(d,elmts(1,4,id,n),FM)
        p = filter_filter(d,elmts(1,5,id,n),FM)
        varpi = atan2(h,k)
        node = zero2pi(atan2(p,q))
        peri = zero2pi(varpi-node)
        if (filter_direct) then
          e   = filter_filter(d,elmts(1,6,id,n),FM)
          inc = filter_filter(d,elmts(1,7,id,n),FM)
        else
          e = sqrt(k*k+h*h)
          inc = 2.0d0 * asin(sqrt(q*q+p*p))
        endif

        s_c = filter_filter(d,elmts(1,8,id,n),FM)
        s_s = filter_filter(d,elmts(1,9,id,n),FM)
        sigma = zero2pi(atan2(s_s,s_c))

        if (filter_write) then    ! write output-file?!
          if (filter_write_r8) then
            call io_write_line_r8(iu,id,a,e,inc,node,peri,sigma)
          else
            call io_write_line_r(iu,id,a,e,inc,node,peri,sigma)
          endif
        endif

c  DO store planetary mean elements also (for FMFT filter in io_write_proper.f)
        mean_elmts(1,id) = a
        mean_elmts(2,id) = k
        mean_elmts(3,id) = h
        mean_elmts(4,id) = q
        mean_elmts(5,id) = p
        mean_elmts(6,id) = e          ! I use these for avg and min/max filter in io_write_proper.f
        mean_elmts(7,id) = inc
        mean_elmts(8,id) = sigma      ! I use this for mean resonant elements in io_write_proper.f
      enddo

      do i=1,ntp
        if (istat(i,1).eq.0) then
          id=i
          a = filter_filter(d,elmts(1,1,id,n),FM)
          k = filter_filter(d,elmts(1,2,id,n),FM)
          h = filter_filter(d,elmts(1,3,id,n),FM)
          q = filter_filter(d,elmts(1,4,id,n),FM)
          p = filter_filter(d,elmts(1,5,id,n),FM)
          varpi = atan2(h,k)
          node = zero2pi(atan2(p,q))
          peri = zero2pi(varpi-node)
c  May 12th 2003: calculate e, inc again from h,k,p,q
c  (Because proper elements are calculted from these and
c  we need proper elements to be reproducible!)
          if (filter_direct) then
            e   = filter_filter(d,elmts(1,6,id,n),FM)
            inc = filter_filter(d,elmts(1,7,id,n),FM)
          else
            e = sqrt(k*k+h*h)
            inc = 2.0d0 * asin(sqrt(q*q+p*p))
          endif

          s_c = filter_filter(d,elmts(1,8,id,n),FM)
          s_s = filter_filter(d,elmts(1,9,id,n),FM)
          sigma = zero2pi(atan2(s_s,s_c))

          if (filter_write) then
            if (filter_write_r8) then
              call io_write_line_r8(iu,id,a,e,inc,node,peri,sigma)
            else
              call io_write_line_r(iu,id,a,e,inc,node,peri,sigma)
            endif
          endif

c  store output mean elements in array mean_elmts (and update flag)
          mean_elmts(1,id) = a
          mean_elmts(2,id) = k
          mean_elmts(3,id) = h
          mean_elmts(4,id) = q
          mean_elmts(5,id) = p
          mean_elmts(6,id) = e          ! I use these for avg and min/max filter in io_write_proper.f
          mean_elmts(7,id) = inc
          mean_elmts(8,id) = sigma      ! I use this for mean resonant elements in io_write_proper.f
        endif
      enddo

      upflg_mean = .true.

      return
      end       ! filter_elmts_write

