
      subroutine proper_orbfit(t,tstart,tstop,nbod,ntp,istat,oname,iu)
c
c  Computation of synthetic proper elements as defined by
c  Knezevic & Milani (2000)
c
c  Author: Ondrej Chrenko (chrenko@sirrah.troja.mff.cuni.cz)
c  Cite: Yang et al. 2020 (Astronomy & Astrophysics, Volume 643, id.A38, 9 pp.)
c

      include '../swift.inc'
      include 'filter.inc'
      include 'proper.inc'
      include 'cb_flt.inc'
      include 'cb_meanel.inc'
      include 'cb_propel.inc'

      real*8 t,tstart,tstop
      integer nbod,ntp,iu,istat(NTPMAX,NSTAT)
      character*(*) oname

      integer PN
      save PN

c  temporal variables
      integer i,id,ierr,nbodtmp,j,nfforced
      real*8 tout,capa,e,inc,tmp,freq,period,amp,ph

      real*8 times(PMAXN),k(PMAXN),h(PMAXN)
      real*8 q(PMAXN),p(PMAXN),varpi(PMAXN),capom(PMAXN)
      real*8 fforced(MAX_FMFT_NFREQKNOWN)

c  functions
      real*8 arr_avg

      data PN /1/

c=======================================================================
c  ... executable code

      prop_time(PN) = t

c      do i = 2, nbod
c        id = -i
c        prop_elmts(PN,1,id) = mean_elmts(1,id)
c        prop_elmts(PN,2,id) = mean_elmts(6,id)
c        prop_elmts(PN,3,id) = mean_elmts(7,id)
c      enddo

      do i = 1, ntp
        if (istat(i,1).eq.0) then
          prop_elmts(PN,1,i) = mean_elmts(1,i)
          prop_elmts(PN,2,i) = mean_elmts(2,i)
          prop_elmts(PN,3,i) = mean_elmts(3,i)
          prop_elmts(PN,4,i) = mean_elmts(4,i)
          prop_elmts(PN,5,i) = mean_elmts(5,i)
        endif
      enddo

c  the buffer is filled-up => write the output
      if (t.gt.tstop) then

        tout = t - 0.5d0*prop_win - 0.5d0*FNinibuf*dtfilter

c        call io_open(iu,oname,'append','UNFORMATTED',ierr)

c  2DO only use this until computation for planets is included as well !!!
c  then go back to nbod
        nbodtmp = 0
        if (prop_writer8) then
c          call io_write_hdr_r8(iu,tout,nbod,ntp,istat)
          call io_write_hdr_r8(iu,tout,nbodtmp,ntp,istat)
        else
c          call io_write_hdr_r(iu,tout,nbod,ntp,istat)
          call io_write_hdr_r(iu,tout,nbodtmp,ntp,istat)
        endif

c  2DO implement computation for planets
c        do i = 2, nbod
c          id = -i
c          capa = arr_avg(prop_elmts(1,1,id), PN)
c          e    = arr_avg(prop_elmts(1,2,id), PN)
c          inc  = arr_avg(prop_elmts(1,3,id), PN)
c          if (prop_writer8) then
c            call io_write_line_r8(iu,id,capa,e,inc,0.0d0,0.0d0,0.0d0)
c          else
c            call io_write_line_r(iu,id,capa,e,inc,0.0d0,0.0d0,0.0d0)
c          endif
c        enddo

        do i = 1, ntp
          if (istat(i,1).eq.0) then

c  semimajor axis is simply a running average
            capa = arr_avg(prop_elmts(1,1,i), PN)

            do j=1,PN
              times(j) = (prop_time(j) - 0.5d0*FNinibuf*dtfilter)/365.25d0
c               times(j) = 0.0501206045*1.d6 + (j-1)*1500.d0
            enddo

c  2DO - there is a large hack in OrbFit designed to deal with
c     eccentricities of Astraea. Should we implement it?

            do j=1,PN
              k(j) = prop_elmts(j,2,i)
              h(j) = prop_elmts(j,3,i)
              varpi(j) = 0.d0
            enddo

            nfforced = fmft_nfreqknown_g
            do j=1,nfforced
              fforced(j) = fmft_freqknown(j)
            enddo

            call orbfit_forced(k,h,times,PN,fforced,nfforced)

            call orbfit_argum(k,h,times,PN,varpi,freq)

            period = TWOPI
            call orbfit_prop(varpi,k,h,PN,period,amp,ph)
            e = amp
c            write(*,*) e

            do j=1,PN
              q(j) = prop_elmts(j,4,i)
              p(j) = prop_elmts(j,5,i)
              capom(j) = 0.d0
            enddo

            nfforced = fmft_nfreqknown - fmft_nfreqknown_g
            do j=1,nfforced
              fforced(j) = fmft_freqknown(j+fmft_nfreqknown_g)
            enddo

c            write(*,*) 'Some input values',q(1),q(PN),p(1),p(PN)
c            write(*,*) times(1),times(PN)
c            write(*,*) fforced(1),fforced(nfforced)
c  all 2DOs below have to be applied for eccentricities as well
            call orbfit_forced(q,p,times,PN,fforced,nfforced)

            call orbfit_argum(q,p,times,PN,capom,freq)
c  2DO include an option to print s
            period = TWOPI
            call orbfit_prop(capom,q,p,PN,period,amp,ph)
            inc = 2.d0*asin(amp)
c            write(*,*) inc*DEGRAD
c            stop

            if (prop_writer8) then
              call io_write_line_r8(iu,i,capa,e,inc,0.0d0,0.0d0,0.0d0)
            else
              call io_write_line_r(iu,i,capa,e,inc,0.0d0,0.0d0,0.0d0)
            endif
          endif
        enddo

c        close(iu)

c  shift the prop_elmts buffer back in time
        call proper_shift(tstart,tstop,PN,nbod,ntp,
     :    prop_time,prop_elmts,5)

      endif

      PN = PN + 1

      return
      end


