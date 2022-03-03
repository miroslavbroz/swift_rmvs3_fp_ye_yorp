
      subroutine proper_sigma_2ND(t,tstart,tstop,nbod,ntp,istat,
     :  oname,iu)
c
c  sigma & varpi-varpi_PL \simeq 0 filter (resonant elements)
c  critical angle of a resonance: sigma = (p+q)/q*lambda_PL - p/q*lambda_TP - varpi_TP
c
c  we actually do NOT need the buffer, we must only check the above condition
c  and take THE LAST OSCULATING elements, instead of mean one! (lambda and
c  sigma changes very quickly and we need a sufficient sampling)
c
      include '../swift.inc'
      include 'filter.inc'
      include 'proper_2ND.inc'
      include 'cb_flt_2ND.inc'
      include 'cb_oscel_2ND.inc'
      include 'cb_propel_2ND.inc'

      real*8 t,tstart,tstop
      integer nbod,ntp,iu,istat(NTPMAX,NSTAT)
      character*(*) oname

      integer PN,i1st
      real*8 sigmalast(-PMAXNPL:PMAXNTP)
      save PN,i1st,sigmalast

c  temporal variables
      integer i,j,ierr
      real*8 tout,Omega,varpi,lambda,Omega_PL,varpi_PL,lambda_PL,sigma

c  functions
      integer ang_sgndot
      real*8 zero2pi
      logical ang_eq

      data PN,i1st /1,0/

c=======================================================================
c  ... executable code

      if (i1st.eq.0) then
        do i = -nbod, ntp
          sigmalast(i) = 0.d0
        enddo
        i1st = 1
      endif

      Omega_PL = osc_elmts(4,prop_plid)
      varpi_PL = zero2pi(Omega_PL + osc_elmts(5,prop_plid))
      lambda_PL = zero2pi(varpi_PL + osc_elmts(6,prop_plid))

      do i = -nbod, ntp
        if ((i.le.-2).or.((i.ge.1).and.(istat(i,1).eq.0))) then

          Omega = osc_elmts(4,i)
          varpi = zero2pi(Omega + osc_elmts(5,i))
          lambda = zero2pi(varpi + osc_elmts(6,i))
          sigma = zero2pi(dble(prop_p+prop_q)/prop_q * lambda_PL
     :      - dble(prop_p)/prop_q * lambda
     :      - varpi)

c  if the condition is true, save the last osculating elements
          if (ang_eq(sigma, prop_sigma, prop_dsigma)
     :      .and.((prop_sigmadot.eq.0).or.
     :      (ang_sgndot(sigmalast(i), sigma).eq.prop_sigmadot))
     :      .and.((prop_dvarpi.lt.0.d0).or.
     :      ang_eq(zero2pi(varpi-varpi_PL),prop_varpi,prop_dvarpi))
     :      .and.((prop_dOmega.lt.0.d0).or.
     :      ang_eq(zero2pi(Omega-Omega_PL),prop_Omega,prop_dOmega)))
     :      then

c  save the OSCULATING elements (regardless short-periodic escilations!)
c  we use the 1st position ONLY in prop_elmts buffer
            do j = 1, 6
              prop_elmts(1,j,i) = osc_elmts(j,i)
            enddo

c  debug code
c            write(*,'(i5,1x,7(1x,f22.10),2(1x,i3))')
c     :        i,
c     :        t/365.25d6,
c     :        osc_elmts(1,i),
c     :        osc_elmts(2,i),
c     :        sigma*degrad,
c     :        zero2pi(varpi-varpi_PL)*degrad,
c     :        zero2pi(Omega-Omega_PL)*degrad,
c     :        sigmalast(i)*degrad,
c     :        ang_sgndot(sigmalast(i), sigma),
c     :        prop_sigmadot
          endif

c  save last sigma due to d sigma/d t monitoring
          sigmalast(i) = sigma
        endif
      enddo

c  the buffer is filled-up => write the output
      if (t.gt.tstop) then

        tout = t - 0.5d0*prop_win

        call io_open(iu,oname,'append','UNFORMATTED',ierr)
        if (prop_writer8) then
          call io_write_hdr_r8(iu,tout,nbod,ntp,istat)
        else
          call io_write_hdr_r(iu,tout,nbod,ntp,istat)
        endif

        do i = -nbod, ntp
          if ((i.le.-2).or.((i.ge.1).and.(istat(i,1).eq.0))) then
            if (prop_writer8) then
              call io_write_line_r8(iu,i,
     :          prop_elmts(1,1,i),
     :          prop_elmts(1,2,i),
     :          prop_elmts(1,3,i),
     :          prop_elmts(1,4,i),
     :          prop_elmts(1,5,i),
     :          prop_elmts(1,6,i))
            else
              call io_write_line_r(iu,i,
     :          prop_elmts(1,1,i),
     :          prop_elmts(1,2,i),
     :          prop_elmts(1,3,i),
     :          prop_elmts(1,4,i),
     :          prop_elmts(1,5,i),
     :          prop_elmts(1,6,i))
            endif
          endif
        enddo

        close(iu)

        tstop = tstop + prop_dt
        tstart = tstop - prop_win
      endif

      return
      end


