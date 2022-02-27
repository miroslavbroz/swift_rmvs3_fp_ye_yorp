
      subroutine proper_trojan(t,tstart,tstop,nbod,ntp,istat,
     :  oname,iu)
c
c  a filter suitable for Trojan asteroids: calculates libration frequency
c  first and than the amplitude of semimajor axis and critical angle
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
      integer i,j,id,ierr
      real*8 a_a_(PMAXN), phi_crit(PMAXN), phi360(PMAXN), time(PMAXN)
      real*8 tout, f, d, capD, dummy, phi, lastphi, phiadd, eps, a0, l0

      data PN /1/
      data eps /1.d-8/

      save eps

c=======================================================================
c  ... executable code

      prop_time(PN) = t

      do i = -nbod, ntp
        if ((i.le.-2).or.((i.ge.1).and.(istat(i,1).eq.0))) then
          prop_elmts(PN,1,i) = mean_elmts(1,i)
          prop_elmts(PN,2,i) = mean_elmts(8,i)
        endif
      enddo

c  the buffer is filled-up => write the output
      if (t.ge.(tstop-eps)) then

        tout = t - 0.5d0*prop_win - 0.5d0*FNinibuf*dtfilter

        call io_open(iu,oname,'append','UNFORMATTED',ierr)
        if (prop_writer8) then
          call io_write_hdr_r8(iu,tout,nbod,ntp,istat)
        else
          call io_write_hdr_r(iu,tout,nbod,ntp,istat)
        endif

        do j = 1, PN
          time(j) = prop_time(j)-tstart
        enddo

        do i = -nbod, ntp
          if ((i.le.-2).or.((i.ge.1).and.(istat(i,1).eq.0))) then
            id = i
            if (id.lt.0) id = -(nbod+i+2)                       ! planet no. 2 is the first
            lastphi = 0.d0
            phiadd = 0.d0

c  fit a-a' sequence by LSM and subtract
c  (to compensate an offset of the libration centre for small libration amplitudes)
            do j = 1, PN
              a_a_(j) = prop_elmts(j,1,id) - prop_elmts(j,1,prop_plid)
            enddo

            call lsm(time, a_a_, PN, dummy, a0, dummy)

c  fit also lambda-lambda'-chi sequence
            do j = 1, PN
              phi_crit(j) = prop_elmts(j,2,id) - prop_sigma    ! phi_crit = lambda - lambda' - chi
              if (phi_crit(j).gt.PI) then
                phi_crit(j) = phi_crit(j) - TWOPI
              endif
            enddo

            call lsm(time, phi_crit, PN, dummy, l0, dummy)

c  calculate the sequence of phi values
            do j = 1, PN
              a_a_(j) = a_a_(j)-a0
              phi_crit(j) = phi_crit(j)-l0
              phi = atan2(a_a_(j)/0.2783d0, phi_crit(j))

c  a 360 deg correction due to linear fit
              if ((phi.lt.-0.5d0*PI).and.(lastphi.gt.0.5d0*PI)) then
                phiadd = phiadd + TWOPI
             else if ((phi.gt.0.5d0*PI).and.(lastphi.lt.-0.5d0*PI)) then
                phiadd = phiadd - TWOPI
              endif
              phi360(j) = (phi + phiadd)
              lastphi = phi
            enddo

c  determine the libration frequency by LSM fit
            call lsm(time, phi360, PN, f, dummy, dummy)

c  calculate the amplitude d of a-a' (for the libration frequency) by simple DFT
            call dft(time, a_a_, PN, f, d, dummy)

c  the same for lambda-lambda'
            call dft(time, phi_crit, PN, f, capD, dummy)

c  write the output (with many zeros...)
            f = f*DEGRAD*365.25d0
            if (prop_writer8) then
              call io_write_line_r8(iu,id,d,f,capD,a0,l0,0.0d0)
            else
              call io_write_line_r(iu,id,d,f,capD,a0,l0,0.0d0)
            endif

          endif
        enddo

        close(iu)

c  shift the prop_elmts buffer back in time
        call proper_shift(tstart,tstop,PN,nbod,ntp,
     :    prop_time,prop_elmts,2)

      endif

      PN = PN + 1

      return
      end


