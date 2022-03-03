
      subroutine proper_fmft_2ND(t,tstart,tstop,nbod,ntp,istat,oname,iu)
c
c  FMFT filter, calls subroutine fmft.c (in C) by D. Nesvorny
c
      include '../swift.inc'
      include 'filter.inc'
      include 'proper_2ND.inc'
      include 'cb_flt_2ND.inc'
      include 'cb_meanel_2ND.inc'
      include 'cb_propel_2ND.inc'

      real*8 t,tstart,tstop
      integer nbod,ntp,iu,istat(NTPMAX,NSTAT)
      character*(*) oname

      integer PN,i1st,i2nd
      real*8 data_sep,minfreq,maxfreq
      save PN,data_sep,minfreq,maxfreq

c  temporal variables
      integer i,j,n,id,ierr,fmft_nfreqknown2,nbod2
      real*8 tout,capa,e,inc,peri,node,varpi
      real*8 prop_time_mean,prop_time_sigma
      real*8 input(PMAXN,2),amp,freq,phase
      real*8 g(-PMAXNPL:PMAXNTP), s(-PMAXNPL:PMAXNTP)
      real*8 goutput(-PMAXNPL:PMAXNTP,MAX_FMFT_NFREQ,3)
      real*8 soutput(-PMAXNPL:PMAXNTP,MAX_FMFT_NFREQ,3)
      logical fmft_restart

c  functions
      real*8 zero2pi,arr_avg

      data PN,i1st,i2nd /1,0,0/

c=======================================================================
c  ... executable code

c  1st time through (only a single test every call)
      if (i1st.eq.0) then
        if (i2nd.eq.0) then
          i2nd = 1
        else
c  2nd time through
c  init for FMFT filter (to determine the time-step)
          if (FNdtout.gt.0) then
            data_sep = FNdtout*dtfilter
          else
            data_sep = (t - tstart)
          endif
          data_sep = data_sep/365.25d0

          minfreq = fmft_minf /180.d0/3600.d0*pi * data_sep
          maxfreq = fmft_maxf /180.d0/3600.d0*pi * data_sep

          i1st = 1
        endif
      endif

      if (fmft_alsopl) then
        nbod2 = nbod
      else
        nbod2 = 0
      endif

      prop_time(PN) = t
      do i = 2, nbod2
        id = -i
        prop_elmts(PN,1,id) = mean_elmts(1,id)
        prop_elmts(PN,2,id) = mean_elmts(2,id)
        prop_elmts(PN,3,id) = mean_elmts(3,id)
        prop_elmts(PN,4,id) = mean_elmts(4,id)
        prop_elmts(PN,5,id) = mean_elmts(5,id)
      enddo

      do i = 1, ntp
        if (istat(i, 1).eq.0) then
          prop_elmts(PN,1,i) = mean_elmts(1,i)
          prop_elmts(PN,2,i) = mean_elmts(2,i)
          prop_elmts(PN,3,i) = mean_elmts(3,i)
          prop_elmts(PN,4,i) = mean_elmts(4,i)
          prop_elmts(PN,5,i) = mean_elmts(5,i)
        endif
      enddo

c  the buffer is filled-up (2^N samples!) => write the output
      if (PN.ge.fmft_ndata) then

c-----------------------------------

c  check for the possible restart interruptions (ie. non equidistant data) in the input buffer
        fmft_restart = .false.
        if (fmft_chkrst) then
          prop_time_mean = 0.d0
          do i = 2, PN
            prop_time_mean = prop_time_mean +
     :        (prop_time(i) - prop_time(i-1))
          enddo
          prop_time_mean = prop_time_mean / (PN-1)
          prop_time_sigma = 0.d0
          do i = 2, PN
            prop_time_sigma = prop_time_sigma +
     :        ((prop_time(i) - prop_time(i-1)) - prop_time_mean)**2
          enddo
          prop_time_sigma = sqrt(prop_time_sigma / (PN-2))
c  3-sigma criterion
          do i = 2, PN
            if (abs((prop_time(i) - prop_time(i-1)) - prop_time_mean)
     :        .gt.(3.d0*prop_time_sigma)) then
              fmft_restart = .true.
            endif
          enddo

c  if it is so, drop the current data and wait for the next buffer
          if (fmft_restart) then
            call proper_shift_2ND(tstart,tstop,PN,nbod2,ntp,
     :        prop_time,prop_elmts,5)
            PN = PN + 1

            return      ! STOP HERE

          endif
        endif

        tout = t

        call io_open(iu,oname,'append','UNFORMATTED',ierr)
        if (prop_writer8) then
          call io_write_hdr_r8(iu,tout,nbod2,ntp,istat)
        else
          call io_write_hdr_r(iu,tout,nbod2,ntp,istat)
        endif

c-----------------------------------

c  call FMFT for each PL and for 2 pairs of elements (h and k, p and q)

c  do NOT discard ANY frequencies for planets!
        fmft_nfreqknown2 = fmft_nfreqknown
        fmft_nfreqknown = 0

        do i = 2, nbod2
          id = -i

c  a (semimajor axis) is an usual running average
          capa = arr_avg(prop_elmts(1,1,id), PN)

c  p and q first
          do j = 1, fmft_ndata
            input(j, 1) = prop_elmts(j,2,id)
            input(j, 2) = prop_elmts(j,3,id)
          enddo

          call fmft_call_2ND(input,minfreq,maxfreq,data_sep,
     :      fmft_ming,fmft_maxg,freq,amp,phase)
          e = amp
          varpi = phase
          g(id) = freq
          call fmft_gsout_2ND(goutput, id)

c  k and h second
          do j = 1, fmft_ndata
            input(j, 1) = prop_elmts(j,4,id)
            input(j, 2) = prop_elmts(j,5,id)
          enddo

          call fmft_call_2ND(input,minfreq,maxfreq,data_sep,
     :      fmft_mins,fmft_maxs,freq,amp,phase)
          inc = 2.0d0*asin(amp)
          node = phase
          s(id) = freq
          call fmft_gsout_2ND(soutput, id)

          peri = varpi-node
          if (peri.lt.0.d0) peri = peri + TWOPI

c  write output (if istat.eq.0)
          if (prop_writer8) then
            call io_write_line_r8(iu,id,capa,e,inc,node,peri,0.0d0)
          else
            call io_write_line_r(iu,id,capa,e,inc,node,peri,0.0d0)
          endif

        enddo

        fmft_nfreqknown = fmft_nfreqknown2

c-----------------------------------

c  call FMFT for each TP and for 2 pairs of elements (h and k, p and q)
        do i = 1, ntp
          if (istat(i,1).eq.0) then

c  a (semimajor axis) is an usual running average
            capa = arr_avg(prop_elmts(1,1,i), PN)

c  p and q first
            do j = 1, fmft_ndata
              input(j, 1) = prop_elmts(j,2,i)
              input(j, 2) = prop_elmts(j,3,i)
            enddo

            call fmft_call_2ND(input,minfreq,maxfreq,data_sep,
     :        fmft_ming,fmft_maxg,freq,amp,phase)
            e = amp
            varpi = phase
            g(i) = freq
            call fmft_gsout_2ND(goutput, i)

c  k and h second
            do j = 1, fmft_ndata
              input(j, 1) = prop_elmts(j,4,i)
              input(j, 2) = prop_elmts(j,5,i)
            enddo

            call fmft_call_2ND(input,minfreq,maxfreq,data_sep,
     :        fmft_mins,fmft_maxs,freq,amp,phase)
            inc = 2.0d0*asin(amp)
            node = phase
            s(i) = freq
            call fmft_gsout_2ND(soutput, i)

            peri = varpi-node
            if (peri.lt.0.d0) peri = peri + TWOPI

c  write output (if istat.eq.0)
            if (prop_writer8) then
              call io_write_line_r8(iu,i,capa,e,inc,node,peri,0.0d0)
            else
              call io_write_line_r(iu,i,capa,e,inc,node,peri,0.0d0)
            endif

          endif
        enddo

        close(iu)

c-----------------------------------

c  optionally write frequencies
        if (fmft_writegs) then
          call io_write_fmft_2ND(tout,nbod2,ntp,istat,g,s,iu,'append')
        endif

c  optionally write frequencies, amplitudes and phases
        if (fmft_propgs) then
          call io_write_fmft_2ND_propgs(tout,nbod2,ntp,istat,
     :      goutput,soutput,iu,'append')
        endif

c  shift the prop_elmts buffer back in time
        call proper_shift_2ND(tstart,tstop,PN,nbod2,ntp,
     :    prop_time,prop_elmts,5)

      endif   ! filled buffer

      PN = PN + 1

      return
      end


