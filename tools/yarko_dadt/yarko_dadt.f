c**********************************************************************
c DADT_YARKO.F
c**********************************************************************
c Computes analytic estimates of linear drift rates da/dt [AU/Myr]
c for both Yarkovsky diurnal and seasonal effects. Input files tp.in,
c th.in and spin.in are the same as for swift_rmvsy integrator.
c
c Author: Miroslav Broz, miroslav.broz@email.cz
c Date: Jan 30th 2003
     
      include 'swift.inc'
      include 'yarko.inc'
      include 'spin.inc'

      real*8 xht(NTPMAX),yht(NTPMAX),zht(NTPMAX)
      real*8 vxht(NTPMAX),vyht(NTPMAX),vzht(NTPMAX)

      real*8 mass(NPLMAX),j2rp2,j4rp4
      real*8 xh(NPLMAX),yh(NPLMAX),zh(NPLMAX)
      real*8 vxh(NPLMAX),vyh(NPLMAX),vzh(NPLMAX)

      integer istat(NTPMAX,NSTAT),i1st
      integer nbod,ntp,nleft
      integer iflgchk,iub,iuj,iud,iue,iuf,iur,iudi
      real*8 rstat(NTPMAX,NSTATR)

      real*8 t0,tstop,dt,dtout,dtdump,dtfilter,dtproper,dtreorient,
     :  dtdisrupt,dtyorp,dtyorpout
      real*8 t,tout,tdump,tfrac,eoff,tfilter,tproper,treorient,
     :  tdisrupt,tyorp,tyorpout

      real*8 rmin,rmax,rmaxu,qmin,rplsq(NPLMAX)
      logical*2 lclose 

      character*80 outfile,infile,fopenstat,outfilterfile,outproperfile

c  reorientation and disruption variables
      integer iseed
      logical reoriented, disrupted
      real*8 tau_reor(NTPMAX),tau_disr(NTPMAX),beta_1,
     :  omega_1,omega_2,a_mean1, a_mean2

c  name of this driver
      character*16 DRIVER
      common /drivername/ DRIVER

c  temorary variables
      real*8 da, da_diurnal(NTPMAX), da_seasonal(NTPMAX),
     &   elmts(6,NTPMAX), irht(NTPMAX), ir3ht(NTPMAX),
     &   obliq(NTPMAX), gm, obliquity
      integer ialpha, i

c  main ----------------------------------------------------------

c  print versiion
      write(*,*) '----------------------------------------------------'
      write(*,*) 'YARKO_DADT version Jul 8th 2010'
      write(*,*) 'Author: Miroslav Broz, Charles University'
      write(*,*) '        e-mail: miroslav.broz@email.cz'
      write(*,*) '        http://sirrah.troja.mff.cuni.cz/~mira/mp/'
      write(*,*) '----------------------------------------------------'

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

c Read spin axes data file
      write(*,*) 'Enter name of spin axes data file : '
      read(*,999) infile
      call io_init_spin(infile,ntp,iseed)

c Read Yarkovsky parameters from file
      write(*,*) 'Enter name of Yarkovsky data file : '
      read(*,999) infile
      call io_init_yarko(infile,ntp)

c Read collisional parameters from file
      write(*,*) 'Enter name of collision data file : '
      read(*,999) infile
      call io_init_collision(infile,nbod,ntp,
     &  dtreorient,dtdisrupt,tau_reor,tau_disr,beta_1,
     &  omega_1,omega_2,a_mean1,a_mean2,iflgchk)

c  calculate osculating orbital elements
      gm = mass(1)
      do i = 1, ntp
        call orbel_xv2el(xht(i),yht(i),zht(i),vxht(i),vyht(i),vzht(i),
     &    gm,ialpha,elmts(1,i),elmts(2,i),elmts(3,i),elmts(4,i),
     &    elmts(5,i),elmts(6,i))
      enddo

      call getacch_ir3(ntp,1,xht,yht,zht,ir3ht,irht)
 
c  calculate obliquities of spin axes

      do i = 1, ntp
        obliq(i) = obliquity(xht(i),yht(i),zht(i),vxht(i),vyht(i),
     :    vzht(i),s(1,i),s(2,i),s(3,i))
      enddo

c  calculate parameters for getacc_yarko subroutine (Yarkovsky diurnal/seasonal force)
      call yarko_omega(ntp,istat)
c      call yarko_seasonal(ntp,mass,istat)

c  calculate drift rates da/dt [AU/Myr]
      call diurnal_dadt(nbod,ntp,mass,xh,yh,zh,xht,yht,zht,istat,irht,
     :  elmts,obliq,da_diurnal)
      call seasonal_dadt(ntp,mass,istat,elmts,da_seasonal)

c  print output
      write(*,10)
10    format("a [AU] e [] inc [deg]  R [m]      obliq. [deg]  T [h]",
     :  "        diurnal  seasonal  da/dt [AU/Myr] ",/,
     :  "========================================================",
     :  "======================================")

      do i = 1, ntp
        da = da_diurnal(i) + da_seasonal(i)
        write(*,20) elmts(1,i), elmts(2,i), elmts(3,i)*DEGRAD,
     :     R(i), obliq(i)*DEGRAD, TWOPI/omega(i)/3600.d0,
     :     da_diurnal(i), da_seasonal(i), da
20      format(f6.3,1x,f5.3,1x,f6.2,2x,f14.6,1x,f7.2,2x,f13.8,
     :    3(1x,f11.8))
      enddo

      end
c---------------------------------------------------------------------

