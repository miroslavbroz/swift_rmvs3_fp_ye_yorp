c**********************************************************************
c SWIFT_MVS2_FP.F
c**********************************************************************
c Integrator of the 2nd order SBAB2 (see Laskar & Robutel, 2000).
c NO close encounters. To run, need 4 input files.
c The code prompts for the file names, but examples are:
c   param.in, pl.in, tp.in, filter.in, proper.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    8/25/94
c Last revision: 12/27/96

c Remarks: Added filter/decimation process, real*8 output.
c   on-line calculation of proper elements, 2nd filter.
c Modified by: Miroslav Broz, miroslav.broz@email.cz
c Date: Aug 19th 2008
     
      include 'swift.inc'
      include '../yarko/yarko.inc'

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 mass(NPLMAX),j2rp2,j4rp4
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

      integer istat(NTPMAX,NSTAT),i1st
      integer nbod,ntp,nleft
      integer iflgchk,iub,iuj,iud,iue,iuf,iup,iuf2,iup2
      real*8 rstat(NTPMAX,NSTATR)

      real*8 t0,tstop,dt,dtout,dtdump,dtfilter,dtproper,
     &  eps
      real*8 t,tout,tdump,tfrac,eoff,tfilter,tproper

      real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
      logical lclose 

      character*80 outfile,infile,fopenstat,outfilterfile,outproperfile

      character*32 DRIVER
      common /drivername/ DRIVER

c  main ----------------------------------------------------------
c...    print version number
      DRIVER = "swift_mvs2_fp"
      call util_version

      use_yarko = .false.
      is_forward = .true.

c-----------------------------------------------------------------------
c
c Read 7 input files
c

c Get data for the run and the test particles
      write(*,*) 'Enter name of parameter data file : '
      read(*,999) infile
999   format(a)
      call io_init_param(infile,t0,tstop,dt,dtout,dtdump,
     &  iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
      write(*,*) ' '
      write(*,*) 'Enter name of planet data file : '
      read(*,999) infile
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
     &  iflgchk)

c Prompt the name and read the proper parameters file
      write(*,*) 'Enter name of proper parameters file : '
      read(*,999) infile
      call io_init_proper(infile,dtproper,outproperfile,
     :  iflgchk)

c Initialize initial time and times for first output and first dump
      t = t0
      tout = t0 + dtout
      tdump = t0 + dtdump
      tfilter = t0 + dtfilter
      tproper = t0 + dtproper
      eps = 1.e-8

c io unit numbers
      iub = 20
      iuj = 30
      iud = 40
      iue = 60
      iuf = 70
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
  
      nleft = ntp
      i1st = 0

c=======================================================================

c
c  Here's the BIG integrator loop
c
      write(*,*) ' ************** MAIN LOOP ****************** '

      do while ((t.le.tstop-eps).and.((ntp.eq.0).or.(nleft.gt.0)))

        call step_kdk2(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &    xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &    vzht,istat,rstat,dt)

        t = t + dt

c-----------------------------------------------------------------------
c
c  Discard criteria
c
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
        if (t.ge.tout-eps) then 

          if (btest(iflgchk,8)) then  ! bit 8 is set (real*8 dump)
            call io_write_frame_r8(t,nbod,ntp,mass,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &        outfile,iub,fopenstat)
          else if (btest(iflgchk,0)) then  ! bit 0 is set
            call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &        iub,fopenstat)
          else if (btest(iflgchk,1)) then  ! bit 1 is set
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
        if (t.ge.tdump-eps) then

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

        if (t.ge.tfilter-eps) then

          call io_write_filter(t,nbod,ntp,mass,xh,yh,zh,
     &      vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &      outfilterfile,iuf,fopenstat)

          tfilter = tfilter + dtfilter
        endif

c-----------------------------------------------------------------------
c
c Proper elements
c

c Check, if osculating or mean elemens have changed, fill the proper elements buffer and eventually write the output.

        if (t.ge.tproper-eps) then

          call io_write_proper(t,nbod,ntp,istat,outproperfile,iup,
     :      fopenstat)

          tproper = tproper + dtproper
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

      end    ! swift_mvs2_fp2.f
c---------------------------------------------------------------------


