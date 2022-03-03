c**********************************************************************
c SWIFT_MVS2_FP2.F
c**********************************************************************
c Integrator of the 2nd order SBAB2 (see Laskar & Robutel, 2000).
c NO close encounters. To run, need 4 input files.
c The code prompts for the file names, but examples are: param.in,
c pl.in, tp.in, filter.in, proper.in
c
c Authors:  Hal Levison \& Martin Duncan
c Date:    8/25/94
c Last revision: 12/27/96
c
c Remarks: Added filter/decimation process, real*8 output.
c   on-line calculation of proper elements, 2nd filter.
c Modified by: Miroslav Broz, miroslav.broz@email.cz
c Date: Aug 19th 2008
     
      include 'swift.inc'

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 mass(NPLMAX),j2rp2,j4rp4
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

      integer istat(NTPMAX,NSTAT),i1st
      integer nbod,ntp,nleft
      integer iflgchk,iub,iuj,iud,iue,iuf,iur,iui,iup,iuf2,iup2
      real*8 rstat(NTPMAX,NSTATR)

      real*8 t0,tstop,dt,dtout,dtdump,dtfilter,dtproper,
     &  dtfilter2,dtproper2,eps
      real*8 t,tout,tdump,tfrac,eoff,tfilter,tproper,
     &  tfilter2,tproper2

      real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
      logical lclose 

      character*80 outfile,inparfile,inplfile,intpfile,fopenstat,
     &  infilterfile,outfilterfile,inproperfile,outproperfile,
     &  infilterfile2,outfilterfile2,inproperfile2,outproperfile2

      character*32 DRIVER
      common /drivername/ DRIVER

c  main ----------------------------------------------------------
c...    print version number
      DRIVER = "swift_mvs2_fp2"
      call util_version

c Get data for the run and the test particles
      write(*,*) 'Enter name of parameter data file : '
      read(*,999) inparfile
999   format(a)
      call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
     &  iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c Prompt and read name of planet data file
      write(*,*) ' '
      write(*,*) 'Enter name of planet data file : '
      read(*,999) inplfile
      call io_init_pl(inplfile,lclose,iflgchk,nbod,mass,xh,yh,zh,
     &  vxh,vyh,vzh,rplsq,j2rp2,j4rp4)

c Get data for the run and the test particles
      write(*,*) 'Enter name of test particle data file : '
      read(*,999) intpfile
      call io_init_tp(intpfile,ntp,xht,yht,zht,vxht,vyht,
     &  vzht,istat,rstat)

c Prompt the name and read the filter parameters file
      write(*,*) 'Enter name of filter parameters file : '
      read(*,999) infilterfile
      call io_init_filter(infilterfile,nbod,ntp,dtfilter,
     :  outfilterfile,iflgchk)

c Prompt the name and read the proper parameters file
      write(*,*) 'Enter name of proper parameters file : '
      read(*,999) inproperfile
      call io_init_proper(inproperfile,dtproper,outproperfile,
     :  iflgchk)

c Prompt the name and read the 2nd filter parameters file
      write(*,*) 'Enter name of filter parameters file : '
      read(*,999) infilterfile2
      call io_init_filter_2ND(infilterfile2,nbod,ntp,dtfilter2,
     :  outfilterfile2,iflgchk)

c Prompt the name and read the 2nd proper parameters file
      write(*,*) 'Enter name of proper parameters file : '
      read(*,999) inproperfile2
      call io_init_proper_2ND(inproperfile2,dtproper2,outproperfile2,
     :  iflgchk)

c Initialize initial time and times for first output and first dump
      t = t0
      tout = t0 + dtout
      tdump = t0 + dtdump
      tfilter = t0 + dtfilter
      tproper = t0 + dtproper
      tfilter2 = t0 + dtfilter2
      tproper2 = t0 + dtproper2
      eps = 1.e-8

c io unit numbers
      iub = 20
      iuj = 30
      iud = 40
      iue = 60
      iuf = 70
      iur = 80
      iui = 90
      iup = 75
      iuf2 = 76
      iup2 = 77

c  initial io write ----------------------------------------------
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

c***************here`s the big loop *************************************
      write(*,*) ' ************** MAIN LOOP ****************** '

      do while ((t.le.tstop-eps).and.((ntp.eq.0).or.(nleft.gt.0)))

        call step_kdk2(i1st,t,nbod,ntp,mass,j2rp2,j4rp4,
     &    xh,yh,zh,vxh,vyh,vzh,xht,yht,zht,vxht,vyht,
     &    vzht,istat,rstat,dt)

        t = t + dt

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

c if it is time, output orb. elements, 
        if(t.ge.tout-eps) then 

          if (btest(iflgchk,8)) then          ! bit 8 is set (real*8 dump)
            call io_write_frame_r8(t,nbod,ntp,mass,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &        outfile,iub,fopenstat)
          else if (btest(iflgchk,0)) then        ! bit 0 is set
            call  io_write_frame(t,nbod,ntp,mass,xh,yh,zh,vxh,
     &        vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,outfile,
     &        iub,fopenstat)
          else if (btest(iflgchk,1)) then        ! bit 1 is set
            call io_write_frame_r(t,nbod,ntp,mass,xh,yh,zh,
     &        vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &        outfile,iub,fopenstat)
          endif

          tout = tout + dtout
        endif

c If it is time, do a dump
        if(t.ge.tdump-eps) then

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

c If it is time, save data to the filter buffer, if buffer is filled up,
c write out filtered orbital elements of test particles.

        if (t.ge.tfilter-eps) then
          call io_write_filter(t,nbod,ntp,mass,xh,yh,zh,
     &      vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &      outfilterfile,iuf,fopenstat)

          tfilter = tfilter + dtfilter
        endif

c Check, if osculating or mean elemens have changed, fill the proper elements
c buffer and eventually write the output.

        if (t.ge.tproper-eps) then
          call io_write_proper(t,nbod,ntp,istat,outproperfile,iup,
     :      fopenstat)

          tproper = tproper + dtproper
        endif

c The same for the 2nd filter.

        if (t.ge.tfilter2-eps) then
          call io_write_filter_2ND(t,nbod,ntp,mass,xh,yh,zh,
     &      vxh,vyh,vzh,xht,yht,zht,vxht,vyht,vzht,istat,
     &      outfilterfile2,iuf2,fopenstat)

          tfilter2 = tfilter2 + dtfilter2
        endif

        if (t.ge.tproper2-eps) then
          call io_write_proper_2ND(t,nbod,ntp,istat,outproperfile2,iup2,
     :      fopenstat)

          tproper2 = tproper2 + dtproper2
        endif

      enddo
c********** end of the big loop from time 't0' to time 'tstop'

c  Do a final dump for possible resumption later 

      call io_dump_pl('dump_pl.dat',nbod,mass,xh,yh,zh,
     &  vxh,vyh,vzh,lclose,iflgchk,rplsq,j2rp2,j4rp4)
      call io_dump_tp('dump_tp.dat',ntp,xht,yht,zht,
     &  vxht,vyht,vzht,istat,rstat)
      call io_dump_param('dump_param.dat',t,tstop,dt,dtout,
     &  dtdump,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile)

      call util_exit(0)
      end    ! swift_mvs2_fp2.f
c---------------------------------------------------------------------

