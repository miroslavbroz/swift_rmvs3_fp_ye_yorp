c**********************************************************************
c SWIFT_RMVS3_FP_YE_YORP.F
c**********************************************************************
c Original integrator RMVS3.
c Includes close encounters. To run, need 9 input files.
c The code prompts for the file names, but examples are:
c   param.in, pl.in, tp.in, filter.in, proper.in, spin.in, yarko.in,
c   yorp.in, collision.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    8/25/94
c Last revision: 12/27/96

c Remarks: Added thermal effects, YORP, filter/decimation process,
c   real*8 output and on-line calculation of proper elements.
c Modified by: Miroslav Broz, miroslav.broz@email.cz
c Date: Jan 19th 2010
     
      include 'swift.inc'

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 mass(NPLMAX),j2rp2,j4rp4
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

      integer istat(NTPMAX,NSTAT),i1st
      integer nbod,ntp,nleft
      integer iflgchk,iub,iuj,iud,iue,iuf,iur,iui,iup
      real*8 rstat(NTPMAX,NSTATR)

      real*8 t0,tstop,dt,dtout,dtdump,dtfilter,dtproper,dtreorient,
     :  dtdisrupt,dtyorp,dtyorpout
      real*8 t,tout,tdump,tfrac,eoff,tfilter,tproper,treorient,
     :  tdisrupt,tyorp,tyorpout

      real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
      logical lclose 

      character*80 outfile,infile,fopenstat,outfilterfile,outproperfile

c  reorientation and disruption variables
      integer iseed
      logical reoriented, disrupted
      real*8 tau_reor(NTPMAX),tau_disr(NTPMAX),beta_1,
     :  omega_1,omega_2,a_mean1, a_mean2

c  name of this driver
      character*32 DRIVER
      common /drivername/ DRIVER

cc  OMP timings, debugging
c!$    real *8 omp_get_wtime
c!$    real *8 tbegin, tend
c!$    common /timing/ tbegin, tend

c!$    tbegin = omp_get_wtime()

c...    print version number
      DRIVER = "swift_rmvs3_fp_ye_YORP"
      call util_version

c-----------------------------------------------------------------------
c
c Read 9 input files
c

c Get data for the run and the test particles
      write(*,*) 'Enter name of parameter data file : '
      read(*,999) infile
      call io_init_param(infile,t0,tstop,dt,dtout,dtdump,
     &  iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
      write(*,*) ' '
      write(*,*) 'Enter name of planet data file : '
      read(*,999) infile
999   format(a)
      call io_init_pl(infile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &  vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
      write(*,*) 'Enter name of test particle data file : '
      read(*,999) infile
      call io_init_tp(infile,ntp,xht,yht,zht,vxht,vyht,
     &  vzht,istat,rstat)

c Prompt the name and read the filter parameters file
      write(*,*) 'Enter name of filter parameters file : '
      read(*,999) infile
      call io_init_filter(infile,nbod,ntp,dtfilter,outfilterfile,
     :  iflgchk)

c Prompt the name and read the proper parameters file
      write(*,*) 'Enter name of proper parameters file : '
      read(*,999) infile
      call io_init_proper(infile,dtproper,outproperfile,
     :  iflgchk)

c Read spin axes data file
      write(*,*) 'Enter name of spin axes data file : '
      read(*,999) infile
      call io_init_spin(infile,ntp,iseed)

c Read Yarkovsky parameters from file
      write(*,*) 'Enter name of Yarkovsky data file : '
      read(*,999) infile
      call io_init_yarko(infile,ntp)

c Read YORP parameters from file
      write(*,*) 'Enter name of YORP data file : '
      read(*,999) infile
      call io_init_yorp(infile,ntp,dtyorp,dtyorpout)

c Read collisional parameters from file
      write(*,*) 'Enter name of collision data file : '
      read(*,999) infile
      call io_init_collision(infile,nbod,ntp,
     &  dtreorient,dtdisrupt,tau_reor,tau_disr,beta_1,
     &  omega_1,omega_2,a_mean1,a_mean2,iflgchk)

c Initialize initial time and times for first output and first dump
      t = t0
      tout = t0 + dtout
      tdump = t0 + dtdump
      tfilter = t0 + dtfilter
      tproper = t0 + dtproper
      tdisrupt = t0 + dtdisrupt
      treorient = t0 + dtreorient
      tyorp = t0 + dtyorp
      tyorpout = t0 + dtyorpout

c io unit numbers
      iub = 20
      iuj = 30
      iud = 40
      iue = 60
      iuf = 70
      iur = 80
      iui = 90
      iup = 75

c-----------------------------------------------------------------------
c
c  Initial io write
c
      if(btest(iflgchk,8))  then ! bit 8 is set
        call io_write_frame_r8(t0,nbod,ntp,mass,xh,yh,zh,
     &    vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,
     &    fopenstat)
      else if(btest(iflgchk,0))  then ! bit 0 is set
        call io_write_frame(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &    xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
      else if(btest(iflgchk,1))  then ! bit 1 is set
        call io_write_frame_r(t0,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &    xht,yht,zht,vxht,vyht,vzht,istat,outfile,iub,fopenstat)
      endif
      if(btest(iflgchk,2))  then    ! bit 2 is set
        eoff = 0.0d0
        call anal_energy_write(t0,nbod,mass,j2rp2,j4rp4,xh,yh,zh,vxh,
     &    vyh,vzh,iue,fopenstat,eoff)
      endif
      if(btest(iflgchk,3))  then    ! bit 3 is set
        call anal_jacobi_write(t0,nbod,ntp,mass,xh,yh,zh,vxh,
     &    vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,iuj,fopenstat)
      endif

c...  must initize discard io routine
      if(btest(iflgchk,4))  then ! bit 4 is set
        call io_discard_write(0,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &    vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &    'discard.out',fopenstat,nleft)
      endif

c  do initial dump of filtered elements (we need to initialize elmts array
c  for getacc_tp.f subroutine)
      call io_write_filter(t,nbod,ntp,mass,xh,yh,zh,
     &  vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &  outfilterfile,iuf,fopenstat)
  
c  calculate parameters for getacc_yarko subroutine (Yarkovsky diurnal/seasonal force)
      call yarko_omega(ntp,istat)
      call yarko_seasonal(ntp,mass,istat)

c  dump spin axes orientations for easier restart
      call io_dump_spin('dump_spin.dat',ntp,iseed)

      if (btest(iflgchk,6)) then ! bit 6 is set
        call reorient_write('reorient.out',t,ntp,istat,iur,fopenstat)
      endif

      nleft = ntp
      i1st = 0

c=======================================================================

c
c  Here's the BIG integrator loop
c
      write(*,*) ' ************** MAIN LOOP ****************** '

      do while ((t.le.tstop).and.((ntp.eq.0).or.(nleft.gt.0)))

        call rmvs3_step(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &    xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &    vzht,istat,rstat,dt)

        t = t + dt

c-----------------------------------------------------------------------
c
c  Disruption/discard criteria
c
        if (t.ge.tdisrupt) then	
          call disrupt(t,ntp,dtdisrupt,tau_disr,iseed,istat,rstat,
     :      disrupted)
          if (btest(iflgchk,7).and.disrupted) then	! bit 7 is set
            call disrupt_write('disrupt.out',t,ntp,istat,rstat,iui,
     :        fopenstat)
          endif

          tdisrupt = tdisrupt + dtdisrupt
        endif

        call discard_meana(t,ntp,a_mean1,a_mean2,istat,rstat)

        if(btest(iflgchk,4))  then    ! bit 4 is set
          call discard(t,dt,nbod,ntp,mass,xh,yh,zh,vxh,vyh,vzh,
     &      xht,yht,zht,vxht,vyht,vzht,rmin,rmax,rmaxu,
     &      qmin,lclose,rplsq,istat,rstat)
          call io_discard_write(1,t,nbod,ntp,xh,yh,zh,vxh,vyh,
     &      vzh,xht,yht,zht,vxht,vyht,vzht,istat,rstat,iud,
     &      'discard.out',fopenstat,nleft)
        else
          nleft = ntp
        endif

c-----------------------------------------------------------------------
c
c Output osculating elements
c
        if (t.ge.tout) then 

          if (btest(iflgchk,8)) then  	! bit 8 is set (real*8 dump)
            call io_write_frame_r8(t,nbod,ntp,mass,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &        outfile,iub,fopenstat)
          else if (btest(iflgchk,0)) then	! bit 0 is set
            call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &        iub,fopenstat)
          else if (btest(iflgchk,1)) then	! bit 1 is set
            call io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &        outfile,iub,fopenstat)
          endif

          tout = tout + dtout
        endif

c-----------------------------------------------------------------------
c
c Dump osculating elements
c
        if (t.ge.tdump) then

          tfrac = (t-t0)/(tstop-t0)
          write(*,998) t,tfrac,nleft
 998      format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &      ': Number of active tp =',i4)
          call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     &      vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
          call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &      vxht,vyht,vzht,istat,rstat)
          call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &      dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)
          tdump = tdump + dtdump

          if(btest(iflgchk,2))  then    ! bit 2 is set
            call anal_energy_write(t,nbod,mass,j2rp2,j4rp4,
     &        xh,yh,zh,vxh,vyh,vzh,iue,fopenstat,eoff)
          endif
          if(btest(iflgchk,3))  then    ! bit 3 is set
            call anal_jacobi_write(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,2,
     &        iuj,fopenstat)
          endif

        endif

c-----------------------------------------------------------------------
c
c Mean elements
c

c Save data to the filter buffer, if buffer is filled up, write mean orbital elements.
c Recalculate parameters for seasonal Yarkovky effect, which depend on mean elements!

        if (t.ge.tfilter) then

          call io_write_filter(t,nbod,ntp,mass,xh,yh,zh,
     &      vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &      outfilterfile,iuf,fopenstat)

          call yarko_seasonal(ntp,mass,istat)

          tfilter = tfilter + dtfilter
        endif

c-----------------------------------------------------------------------
c
c Proper elements
c

c Check, if osculating or mean elemens have changed, fill the proper elements buffer and eventually write the output.

        if (t.ge.tproper) then

          call io_write_proper(t,nbod,ntp,istat,outproperfile,iup,
     :      fopenstat)

          tproper = tproper + dtproper
        endif

c-----------------------------------------------------------------------
c
c Evolve spin axes due to YORP effect
c

c Recompute all Yarkovsky parameters depending on omega!

        if (t.ge.tyorp) then
          call yorp_evolve(t,ntp,dtyorp,mass,xht,yht,zht,vxht,vyht,vzht,
     :      istat)
          call omega_crit(ntp,omega_1,omega_2,iseed,istat)
          call yarko_omega(ntp,istat)

          call io_dump_spin('dump_spin.dat',ntp,iseed)

          tyorp = tyorp + dtyorp
        endif

        if (t.ge.tyorpout) then
          if (btest(iflgchk,6)) then	! bit 6 is set
            call reorient_write('reorient.out',t,ntp,istat,iur,
     :        fopenstat)
          endif

          tyorpout = tyorpout + dtyorpout
        endif

c-----------------------------------------------------------------------
c
c Reorient spin axes due to collisions
c
        if (t.ge.treorient) then
          call reorient(t,ntp,dtreorient,tau_reor,beta_1,
     :      omega_1,omega_2,iseed,istat,reoriented)

          if ((btest(iflgchk,6)).and.(reoriented)) then	! bit 6 is set
            call reorient_write('reorient.out',t,ntp,istat,iur,
     :        fopenstat)
          endif

          if (reoriented) then
            call io_dump_spin('dump_spin.dat',ntp,iseed)
          endif

          treorient = treorient + dtreorient
        endif

      enddo

c end of the big loop from time 't0' to time 'tstop'

c=======================================================================

c  Do a final dump for possible resumption later 

      call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     &  vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
      call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &  vxht,vyht,vzht,istat,rstat)
      call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &  dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

      call io_close(iub,iuf,iup)
      call util_exit(0)

      end    ! swift_mvs2_fp_ye.f
c---------------------------------------------------------------------


