
      subroutine proper_sigma2(t,tstart,tstop,nbod,ntp,istat,
     :  oname,iu)
c
c  mean sigma & mean varpi-varpi_PL \simeq 0 filter (mean resonant elements)
c  critical angle of a resonance: sigma = (p+q)/q*lambda_PL - p/q*lambda_TP - varpi_TP
c
c  we actually do NOT need the buffer, we must only check the above condition
c  and take THE AVERAGE OF MEAN elements, instead of an osculating one!
c  (lambda and sigma changes very quickly and we need a sufficient sampling)
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

      integer PN,i1st,n(-PMAXNPL:PMAXNTP)
      real*8 sigmalast(-PMAXNPL:PMAXNTP),
     :  alast(-PMAXNPL:PMAXNTP),
     :  elast(-PMAXNPL:PMAXNTP),
     :  ilast(-PMAXNPL:PMAXNTP)
      save PN,i1st,n,sigmalast,alast,elast,ilast

c  temporal variables
      integer i,j,ierr
      real*8 tout,capa,e,inc,h,k,q,p
      real*8 Omega,varpi,lambda,Omega_PL,varpi_PL,lambda_PL,sigma

c  functions
      integer ang_sgndot
      real*8 zero2pi
      logical ang_eq

      data PN,i1st /1,0/

c=======================================================================
c  ... executable code

      if (i1st.eq.0) then
        do i = -nbod, ntp
          n(i) = 0
          prop_elmts(1,1,i) = 0.d0
          prop_elmts(1,2,i) = 0.d0 
          prop_elmts(1,3,i) = 0.d0
          sigmalast(i) = 0.d0
          alast(i) = 0.d0
          elast(i) = 0.d0
          ilast(i) = 0.d0
        enddo
        i1st = 1
      endif

      k = mean_elmts(2,prop_plid)
      h = mean_elmts(3,prop_plid)
      q = mean_elmts(4,prop_plid)
      p = mean_elmts(5,prop_plid)
      varpi_PL = zero2pi(atan2(h,k))
      Omega_PL = zero2pi(atan2(p,q))

c  store both PLs and TPs (in case of SyMBA, we have PLs only)
      do i = -nbod, ntp
        if ((i.le.-2).or.((i.ge.1).and.(istat(i,1).eq.0))) then

          k = mean_elmts(2,i)
          h = mean_elmts(3,i)
          q = mean_elmts(4,i)
          p = mean_elmts(5,i)
          varpi = zero2pi(atan2(h,k))
          Omega = zero2pi(atan2(p,q))
          sigma = mean_elmts(8,i)       ! this was calculated in io_write_filter subroutine

c  if this condition is true, we are in the representative surface, luckily
          if (ang_eq(sigma, prop_sigma, prop_dsigma)
     :      .and.((prop_sigmadot.eq.0).or.
     :      (ang_sgndot(sigmalast(i), sigma).eq.prop_sigmadot))
     :      .and.((prop_dvarpi.lt.0.d0).or.
     :      ang_eq(zero2pi(varpi-varpi_PL),prop_varpi,prop_dvarpi))
     :      .and.((prop_dOmega.lt.0.d0).or.
     :      ang_eq(zero2pi(Omega-Omega_PL),prop_Omega,prop_dOmega)))
     :      then
            n(i) = n(i) + 1
            prop_elmts(1,1,i) = prop_elmts(1,1,i) + mean_elmts(1,i)
            prop_elmts(1,2,i) = prop_elmts(1,2,i) + sqrt(k*k+h*h)
            prop_elmts(1,3,i) = prop_elmts(1,3,i)
     :        + 2.0d0 * asin(sqrt(q*q+p*p))
          endif

c  save the last sigma due to d sigma/d t monitoring
          sigmalast(i) = sigma
        endif
      enddo

c  the buffer is filled-up => test the resonance conditions
c  and then calculate average and write the output
      if (t.gt.tstop) then

        tout = t - 0.5d0*prop_win - 0.5d0*FNinibuf*dtfilter

        call io_open(iu,oname,'append','UNFORMATTED',ierr)
        if (prop_writer8) then
          call io_write_hdr_r8(iu,tout,nbod,ntp,istat)
        else
          call io_write_hdr_r(iu,tout,nbod,ntp,istat)
        endif

        do i = -nbod, ntp
          if ((i.le.-2).or.((i.ge.1).and.(istat(i,1).eq.0))) then

c  store elements and use them in future timestep, if the conditions is not fulfilled anymore
            if (n(i).gt.0) then
              alast(i) = prop_elmts(1,1,i) / n(i)
              elast(i) = prop_elmts(1,2,i) / n(i)
              ilast(i) = prop_elmts(1,3,i) / n(i)
            endif

            if (prop_writer8) then
              call io_write_line_r8(iu,i,alast(i),elast(i),ilast(i),
     :          0.d0,0.d0,0.0d0)
            else
              call io_write_line_r(iu,i,alast(i),elast(i),ilast(i),
     :          0.d0,0.d0,0.0d0)
            endif

            n(i) = 0               ! THESE WERE MISSING! corrected on Apr 30th 2008
            prop_elmts(1,1,i) = 0.
            prop_elmts(1,2,i) = 0.
            prop_elmts(1,3,i) = 0.

          endif
        enddo

        close(iu)

        tstop = tstop + prop_dt
        tstart = tstop - prop_win

c  we do NOT shift sigmalast(i) => it works all right ONLY for twin=dt

      endif   ! t.gt.tstop

      return
      end


