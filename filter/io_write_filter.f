c
c             Input:
c                 time          ==>  current time (real scalar)
c                 nbod          ==>  number of massive bodies (int scalar)
c                 ntp           ==>  number of test particles (int scalar)
c                 mass          ==>  mass of bodies (real array)
c                 xh,yh,zh      ==>  current position in helio coord 
c                                    (real arrays)
c                 vxh,vyh,vzh   ==>  current velocity in helio coord 
c                                    (real arrays)
c                 xht,yht,zht    ==>  current part position in helio coord 
c                                      (real arrays)
c                 vxht,vyht,vzht ==>  current velocity in helio coord 
c                                        (real arrays)
c                 istat           ==>  status of the test paricles
c                                      (2d integer array)
c                                      istat(i,1) = 0 ==> active:  = 1 not
c                                      istat(i,2) = -1 ==> Danby did not work
c                 oname           ==> output file name (character string) 
c                 iu              ==> unit number to write to
c                 fopenstat       ==>  The status flag for the open 
c                                      statements of the output files.  
c                                          (character*80)
c
c Other filtration parameters in common block /filter/!
c
c Remarks: Based on io_write_frame
c Authors: Miroslav Broz, miroslav.broz@email.cz
c Date:    Dec 5th 1998
c Revision: Dec 8th 2003
c - fmft filter now correctly handles shift of p, q buffers
c Revision: May 12th 2003
c - output filtered elements e and inc are calculated from h,k,p,q
c (in accordance with proper elements filter, which uses h,k,p,q variables)
c Revision: Mar 15th 2002
c - now the filtering is done for a,h=e*sin(peri+node),k,p,q,e,inc
c - corrected startup skip of buffers due to alignment of output data
c - added OpenMP directives for parallelizing compilers
c Revision: Dec 2nd 2002
c - corrected output time (shifted by half of the initial buffer fillup)
c - flt_cb.inc include file with /filter/ common block
c Revision: Mar 18th 2003
c - write of mean elements might be switched-off
c Revision: Apr 8th 2003
c - output could be in real*8 precision
c Revision: Jul 25th 2006
c - added filtering of a slow variable sigma
c Revision: Aug 4th 2008
c - sigma was multiplied by q (due to Trojans)

      subroutine io_write_filter(time,nbod,ntp,mass,xh,yh,zh,vxh,
     &           vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,oname,
     &           iu,fopenstat)

      include '../swift.inc'
      include 'filter.inc'
      include 'cb_flt.inc'
      include 'cb_oscel.inc'

c...  Inputs: 
      integer nbod,ntp,iu
      real*8 mass(nbod),time
      integer istat(NTPMAX,NSTAT)
      real*8 xh(nbod),yh(nbod),zh(nbod)
      real*8 vxh(nbod),vyh(nbod),vzh(nbod)
      real*8 xht(ntp),yht(ntp),zht(ntp)
      real*8 vxht(ntp),vyht(ntp),vzht(ntp)
      character*80 oname,fopenstat
      integer iflgchk

c...  Internals
      integer i,n,id,k
      integer ialpha,ierr,FNskip,FNcnt

      real*8 gm,capa,e,inc,node,peri,capm,varpi,sininc2,
     :  lambda,lambda_PL,sigma,tout
      integer i1st,i2nd    ! =0 first time through; =1 after
      integer FN(FMAXC)         ! index in elmts() buffer, (1;FBUFN)
      save i1st,i2nd,FN,FNskip,FNcnt
c
c  orbital elements are stored in common block /elements/, will be
c  used for computation of seasonal Yarkovsky force in getacc_tp.f;
c  dimensions of the array: i) width of the buffer, ii) 7 orbital
c  elements, iii) PLs and TPS, iv) number of buffer levels
c 
      real*8 elmts(FMAXN,FMAXE,-FMAXNPL:FMAXNTP,FMAXC)
      save elmts        ! VERY IMPORTANT (otherwise omf77 crashes on SIGSEGV)

c  g77/RH60 compiler does NOT include this equivalence => assign idx=FN(1)
c  just before return from this subroutine
c      equivalence (idx,FN(1))
c  another solution -> pass only a small copy of elmts() array (see oscel_cb.inc)

      data i1st,i2nd,FNcnt/0,0,0/

c  function
      real*8 zero2pi

c----
c...  Executable code 

c  first time through
      if (i1st.eq.0) then
c  do this block only once
        if (i2nd.eq.0) then
          if (filter_write) then
            call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
            if(ierr.ne.0) then
              write(*,*) ' SWIFT ERROR: in io_write_filter: '
              write(*,*) '     Could not open binary output file:'
              call util_exit(1)
            endif
            close(iu)
          endif

c  index for filter_N levels of the buffer
          do i = 1, filter_N
            FN(i) = 0
          enddo
c
c  Skip first data due to time-alignment in filtered output.
c
c  FNskip ... how many points should we skip?
          FNskip = (FNinibuf/FNdtout + 1) * FNdtout - FNinibuf
          i2nd = 1
        endif

c  wait until we can skip back to the beginning of the buffer
c  (FNcnt is the counter)
        if (FNcnt.lt.FNskip) then
          FNcnt = FNcnt + 1
          FN(1) = 0
        else
          i1st = 1
        endif
      endif

c  calculate orbital elements for planets and test particles
      t0 = time                 ! remember the epoch of mean anomaly
      upflg_osc = .true.        ! update flag for io_write_proper.f subroutine
      FN(1)=FN(1)+1
      gm=mass(1)
c      write(*,*) 'io_write_filter: time = ', time, ' FN = ', FN(1),
c     :  FN(2), FN(3), FN(4)

!$OMP PARALLEL DO
!$OMP& PRIVATE (i, gm, id, ialpha, capa, e, inc, node, peri, capm,
!$OMP&   varpi, sininc2)
!$OMP& SHARED (nbod, xh, yh, zh, vxh, vyh, vzh, elmts)
      do i = 2, nbod
        gm=mass(1)+mass(i)
        id=-i
        call orbel_xv2el(xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i),
     &    gm,ialpha,capa,e,inc,node,peri,capm)

        elmts(FN(1),1,id,1) = capa
        varpi = peri + node
        elmts(FN(1),2,id,1) = e * cos(varpi)
        elmts(FN(1),3,id,1) = e * sin(varpi)
        sininc2 = sin(0.5d0*inc)
        elmts(FN(1),4,id,1) = sininc2 * cos(node)
        elmts(FN(1),5,id,1) = sininc2 * sin(node)

        elmts(FN(1),6,id,1) = e
        elmts(FN(1),7,id,1) = inc

c  I do need sigma of PLs in case of SyMBA
c  a slow variable sigma = (p+q)/q*lambda_PL - p/q*lambda_TP - varpi_TP
        lambda = zero2pi(varpi + capm)
        lambda_PL = zero2pi(osc_elmts(4,filter_plid)
     :    + osc_elmts(5,filter_plid) + osc_elmts(6,filter_plid))
c        sigma = zero2pi(dble(filter_p+filter_q)/filter_q * lambda_PL
c     :    - dble(filter_p)/filter_q * lambda
c     :    - varpi)
        sigma = zero2pi((filter_p+filter_q) * lambda_PL
     :    - filter_p * lambda - filter_q * varpi)
        elmts(FN(1),8,id,1) = cos(sigma)
        elmts(FN(1),9,id,1) = sin(sigma)

c  save osculating elements
        osc_elmts(1,id) = capa
        osc_elmts(2,id) = e
        osc_elmts(3,id) = inc
        osc_elmts(4,id) = peri
        osc_elmts(5,id) = node
        osc_elmts(6,id) = capm
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
!$OMP& PRIVATE (i, id, ialpha, capa, e, inc, node, peri, capm,
!$OMP&   varpi, sininc2)
!$OMP& SHARED (ntp, xht, yht, zht, vxht, vyht, vzht, gm,
!$OMP&   elmts, osc_elmts)
      do i = 1, ntp
        if(istat(i,1).eq.0) then
          id=i
          call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i),
     &      gm,ialpha,capa,e,inc,node,peri,capm)

          elmts(FN(1),1,id,1) = capa
c  Sep 20th 2000: changed the angle in the definition of h, k from peri to varpi
          varpi = zero2pi(peri + node)
          elmts(FN(1),2,id,1) = e * cos(varpi)
          elmts(FN(1),3,id,1) = e * sin(varpi)
          sininc2 = sin(0.5d0*inc)
          elmts(FN(1),4,id,1) = sininc2 * cos(node)
          elmts(FN(1),5,id,1) = sininc2 * sin(node)
c  Sep 20th 2000: separate filtering of e, inc

          elmts(FN(1),6,id,1) = e
          elmts(FN(1),7,id,1) = inc

c  a slow variable sigma = (p+q)/q*lambda_PL - p/q*lambda_TP - varpi_TP
          lambda = zero2pi(varpi + capm)
          lambda_PL = zero2pi(osc_elmts(4,filter_plid)
     :      + osc_elmts(5,filter_plid) + osc_elmts(6,filter_plid))
          sigma = zero2pi((filter_p+filter_q) * lambda_PL
     :      - filter_p * lambda - filter_q * varpi)
          elmts(FN(1),8,id,1) = cos(sigma)
          elmts(FN(1),9,id,1) = sin(sigma)

c  save osculating elements
          osc_elmts(1,id) = capa
          osc_elmts(2,id) = e
          osc_elmts(3,id) = inc
          osc_elmts(4,id) = peri
          osc_elmts(5,id) = node
          osc_elmts(6,id) = capm
        endif
      enddo
!$OMP END PARALLEL DO

c
c  test if buffers are filled up
c
      do n=1,filter_N-1
c 2DO: this does NOT work for filter_MA.ne.filter_MB
c        write(*,*) 'n = ',n,' FN = ',FN(n),' FBUFN = ',FBUFN
        if (FN(n).ge.FBUFN) then

          FN(n+1)=FN(n+1)+1

c  filter the array
          if (filter_stat(n:n).eq.'A') then
            call filter_elmts(elmts,nbod,ntp,istat,FN,n,
     :        filter_MA,filter_A)
          else  ! if (filter_stat(n:n).eq.'B') then
            call filter_elmts(elmts,nbod,ntp,istat,FN,n,
     :        filter_MB,filter_B)
          endif

c  decime with factor filter_D
          FN(n)=FN(n)-filter_D(n)

c  shift left
          call filter_shift(elmts,nbod,ntp,n,FBUFN,filter_D)
        endif
      enddo
c
c  test if the last buffer is filled up and write filtered elmts
c
      if (FN(filter_N).ge.FBUFN) then
        tout = time - FNinibuf*dtfilter/2.d0

        if (filter_write) then
          call io_open(iu,oname,'append','UNFORMATTED',ierr)
          if (filter_write_r8) then
            call io_write_hdr_r8(iu,tout,nbod,ntp,istat)
          else
            call io_write_hdr_r(iu,tout,nbod,ntp,istat)
          endif
        endif

        n=filter_N
        if (filter_stat(n:n).eq.'A') then
          call filter_elmts_write(iu,elmts,n,nbod,ntp,istat,
     &      filter_A,filter_MA,filter_direct,filter_write,
     &      filter_write_r8)
        else    ! if (filter_stat(n:n).eq.'B') then
          call filter_elmts_write(iu,elmts,n,nbod,ntp,istat,
     &      filter_B,filter_MB,filter_direct,filter_write,
     &      filter_write_r8)
        endif

c  "decime" with factor filter_D
        FN(n)=FN(n)-filter_D(n)

c  shift left
        call filter_shift(elmts,nbod,ntp,n,FBUFN,filter_D)

        if (filter_write) then
          close(iu)
        endif

      endif

      return
      end       ! io_write_filter
c----------------------------------------------------------------------


