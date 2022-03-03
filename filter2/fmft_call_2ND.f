
c****************************************************************
      subroutine fmft_call_2ND(input, minfreq, maxfreq, data_sep,
     :  min_gs, max_gs, freq, amp, phase)
c****************************************************************
c
c  call FMFT subroutine, discard small or known frequencies
c  and return a proper element
c
      include '../swift.inc'
      include 'filter.inc'
      include 'proper_2ND.inc'

      real*8 input(MAX_FMFT_NDATA, 2)   ! input array (either (h, k) or (p, q))
      real*8 minfreq, maxfreq           ! minimum and maximum frequency (arcsec/yr)
      real*8 data_sep                   ! time separation of data
      real*8 min_gs, max_gs             ! required range of frequencies
      real*8 freq, amp, phase           ! proper frequency, amplitude and phase

c  g and s frequencies, amplitudes and phases from fmft_call_2ND subroutine
      real*8 gsoutput(MAX_FMFT_NFREQ,3)
      common/gs_output/gsoutput

      real*8 output(MAX_FMFT_NFREQ, 3*MAX_FMFT_FLAG)
      integer j, k, err
      logical fmft_print

      integer vectadd_n
      real*8 vectadd_x, vectadd_y, vectadd_f
c  functions
      integer fmft

      err = fmft(output, fmft_nfreq, minfreq, maxfreq,
     :  fmft_flag, input, fmft_ndata)

c  check the error flag, but do NOT stop
      if (err.eq.0) then
        write(*,20)
20      format('Unknown error in fmft C-function.')
c        stop
      endif

c optionally, save ALL frequencies, amplitudes and phases, like in propgs output
      if (fmft_propgs) then
        do j = 1, fmft_nfreq
          gsoutput(j, 1) = output(j+1, 3*fmft_flag-1)
     :      *((180.d0/pi)*3600.d0)/data_sep
          gsoutput(j, 2) = output(j+1, 3*fmft_flag)
          phase = output(j+1, 3*fmft_flag+1)
          if (phase.lt.0) phase = phase + TWOPI
          gsoutput(j, 3) = phase
        enddo
      endif

c  default null output values
      freq = 0.d0
      amp = 0.d0

      j = 2
      fmft_print = .false.
      vectadd_n = 0
      vectadd_x = 0.d0
      vectadd_y = 0.d0
      vectadd_f = 0.d0

      do while ((j.le.fmft_nfreq+1).and.(.not.fmft_print))

        freq = output(j, 3*fmft_flag-1)
     :    *((180.d0/pi)*3600.d0)/data_sep
        fmft_print = .true.

c  discard too small frequency
        if (abs(freq).lt.fmft_freqsmall) then
          fmft_print = .false.

c  discard frequencies out of (g or s) range
        else if ((freq.lt.min_gs).or.(freq.gt.max_gs)) then
          fmft_print = .false.

c  discard well known planetary frequences
        else
          do k = 1, fmft_nfreqknown
c            if (abs(freq - fmft_freqknown(k)).lt.fmft_freqtoler) then   ! the tolerance is given as an absolute value (arcsec/yr)
            if (abs(freq - fmft_freqknown(k)).lt.                       ! the tolerance is given as FRACTION of *f* /from Oct 11th 2010/
     :        fmft_freqtoler*abs(fmft_freqknown(k))) then
              fmft_print = .false.
            endif
          enddo
        endif

        if (fmft_print) then
          amp = output(j, 3*fmft_flag)
          phase = output(j, 3*fmft_flag+1)
          if (phase.lt.0) phase = phase + TWOPI

c  eventually sum (vectorially) all relevant frequiencies
          if (fmft_vectadd) then
            vectadd_n = vectadd_n + 1
            vectadd_x = vectadd_x + amp*cos(phase)
            vectadd_y = vectadd_y + amp*sin(phase)
            vectadd_f = vectadd_f + freq

            fmft_print = .false.
          endif
        endif

        j = j + 1
      enddo

c  calculate average frequency and its amplitude
      if (fmft_vectadd) then
        freq = vectadd_f / vectadd_n
        amp = sqrt(vectadd_x**2 + vectadd_y**2)
        phase = atan2(vectadd_x, vectadd_y)
      endif

      return
      end

